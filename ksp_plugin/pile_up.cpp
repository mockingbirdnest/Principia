
#include "ksp_plugin/pile_up.hpp"

#include <functional>
#include <list>
#include <map>
#include <memory>

#include "base/map_util.hpp"
#include "geometry/identity.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/part.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

using base::check_not_null;
using base::FindOrDie;
using base::make_not_null_unique;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::Identity;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::RigidTransformation;
using geometry::Velocity;
using geometry::Wedge;
using physics::DegreesOfFreedom;
using physics::RigidMotion;
using quantities::AngularMomentum;
using quantities::Time;
using quantities::si::Radian;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

PileUp::PileUp(
    std::list<not_null<Part*>>&& parts,
    Instant const& t,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    std::function<void()> deletion_callback)
    : lock_(make_not_null_unique<absl::Mutex>()),
      parts_(std::move(parts)),
      ephemeris_(ephemeris),
      adaptive_step_parameters_(adaptive_step_parameters),
      fixed_step_parameters_(fixed_step_parameters),
      history_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      deletion_callback_(std::move(deletion_callback)) {
  LOG(INFO) << "Constructing pile up at " << this;
  MechanicalSystem<Barycentric, NonRotatingPileUp> mechanical_system;
  for (not_null<Part*> const part : parts_) {
    mechanical_system.AddRigidBody(
        part->rigid_motion(), part->mass(), part->inertia_tensor());
  }
  auto const barycentre = mechanical_system.centre_of_mass();
  history_->Append(t, barycentre);

  angular_momentum_ = mechanical_system.AngularMomentum();

  RigidMotion<Barycentric, NonRotatingPileUp> const barycentric_to_pile_up =
      mechanical_system.LinearMotion().Inverse();
  for (not_null<Part*> const part : parts_) {
    actual_part_rigid_motion_.emplace(
        part, barycentric_to_pile_up * part->rigid_motion());
  }
  psychohistory_ = history_->NewForkAtLast();

  RecomputeFromParts();
}

PileUp::~PileUp() {
  LOG(INFO) << "Destroying pile up at " << this;
  if (deletion_callback_ != nullptr) {
    deletion_callback_();
  }
}

std::list<not_null<Part*>> const& PileUp::parts() const {
  return parts_;
}

void PileUp::SetPartApparentRigidMotion(
    not_null<Part*> const part,
    RigidMotion<RigidPart, ApparentBubble> const& rigid_motion) {
  auto const [_, inserted] =
      apparent_part_rigid_motion_.emplace(part, rigid_motion);
  CHECK(inserted) << "Duplicate part " << part->ShortDebugString() << " at "
                  << rigid_motion;
}

Status PileUp::DeformAndAdvanceTime(Instant const& t) {
  absl::MutexLock l(lock_.get());
  Status status;
  if (psychohistory_->back().time < t) {
    DeformPileUpIfNeeded();
    status = AdvanceTime(t);
    NudgeParts();
  }
  return status;
}

void PileUp::RecomputeFromParts() {
  mass_ = Mass();
  intrinsic_force_ = Vector<Force, Barycentric>();
  intrinsic_torque_ = Bivector<Torque, NonRotatingPileUp>();
  angular_momentum_change_ = Bivector<AngularMomentum, NonRotatingPileUp>();

  for (not_null<Part*> const part : parts_) {
    mass_ += part->mass();

    intrinsic_force_ += part->intrinsic_force();

    RigidMotion<RigidPart, NonRotatingPileUp> const part_motion =
        FindOrDie(actual_part_rigid_motion_, part);
    DegreesOfFreedom<NonRotatingPileUp> const part_dof =
        part_motion({RigidPart::origin, RigidPart::unmoving});
    intrinsic_torque_ +=
        Wedge(part_dof.position() - NonRotatingPileUp::origin,
              Identity<Barycentric, NonRotatingPileUp>()(
                  part->intrinsic_force())) * Radian +
        Identity<Barycentric, NonRotatingPileUp>()(part->intrinsic_torque());

    AngularVelocity<NonRotatingPileUp> const part_angular_velocity =
        part_motion.Inverse().angular_velocity_of_to_frame();
    InertiaTensor<NonRotatingPileUp> part_inertia_tensor =
        part_motion.orthogonal_map()(part->inertia_tensor());
    // KSP makes the inertia tensor vary proportionally to the mass; this
    // corresponds to the body uniformly changing density.
    angular_momentum_change_ +=
        Wedge(part_dof.position() - NonRotatingPileUp::origin,
              part->mass_change() * part_dof.velocity()) * Radian +
        part->mass_change() / part->mass() *
            (part_inertia_tensor * part_angular_velocity);
  }
}

void PileUp::WriteToMessage(not_null<serialization::PileUp*> message) const {
  for (not_null<Part*> const part : parts_) {
    message->add_part_id(part->part_id());
  }
  history_->WriteToMessage(message->mutable_history(),
                           /*forks=*/{psychohistory_});
  for (auto const& pair : actual_part_rigid_motion_) {
    auto const part = pair.first;
    auto const& rigid_motion = pair.second;
    rigid_motion.WriteToMessage(
        &((*message->mutable_actual_part_rigid_motion())[part->part_id()]));
  }
  for (auto const& pair : apparent_part_rigid_motion_) {
    auto const part = pair.first;
    auto const& rigid_motion = pair.second;
    rigid_motion.WriteToMessage(
        &((*message->mutable_apparent_part_rigid_motion())[part->part_id()]));
  }
  adaptive_step_parameters_.WriteToMessage(
      message->mutable_adaptive_step_parameters());
  fixed_step_parameters_.WriteToMessage(
      message->mutable_fixed_step_parameters());
}

not_null<std::unique_ptr<PileUp>> PileUp::ReadFromMessage(
    serialization::PileUp const& message,
    std::function<not_null<Part*>(PartId)> const& part_id_to_part,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    std::function<void()> deletion_callback) {
  std::list<not_null<Part*>> parts;
  for (auto const part_id : message.part_id()) {
    parts.push_back(part_id_to_part(part_id));
  }

  bool const is_pre_cartan = !message.has_adaptive_step_parameters() ||
                             !message.has_fixed_step_parameters();
  bool const is_pre_cesàro = message.history().children().empty();
  bool const is_pre_frege = message.actual_part_degrees_of_freedom_size() > 0 ||
                            message.apparent_part_degrees_of_freedom_size() > 0;
  std::unique_ptr<PileUp> pile_up;
  if (is_pre_cesàro) {
    if (is_pre_cartan) {
      pile_up = std::unique_ptr<PileUp>(new PileUp(
          std::move(parts),
          DefaultPsychohistoryParameters(),
          DefaultHistoryParameters(),
          DiscreteTrajectory<Barycentric>::ReadFromMessage(message.history(),
                                                           /*forks=*/{}),
          /*psychohistory=*/nullptr,
          ephemeris,
          std::move(deletion_callback)));
    } else {
      pile_up = std::unique_ptr<PileUp>(new PileUp(
          std::move(parts),
          Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
              message.adaptive_step_parameters()),
          Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
              message.fixed_step_parameters()),
          DiscreteTrajectory<Barycentric>::ReadFromMessage(message.history(),
                                                           /*forks=*/{}),
          /*psychohistory=*/nullptr,
          ephemeris,
          std::move(deletion_callback)));
    }
    // Fork a psychohistory for compatibility if there is a non-authoritative
    // point.
    if (pile_up->history_->Size() == 2) {
      Instant const history_begin_time = pile_up->history_->front().time;
      pile_up->psychohistory_ =
          pile_up->history_->NewForkWithCopy(history_begin_time);
      pile_up->history_->ForgetAfter(history_begin_time);
    } else {
      pile_up->psychohistory_ = pile_up->history_->NewForkAtLast();
    }
  } else {
    DiscreteTrajectory<Barycentric>* psychohistory = nullptr;
    not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> history =
        DiscreteTrajectory<Barycentric>::ReadFromMessage(
            message.history(),
            /*forks=*/{&psychohistory});
    pile_up = std::unique_ptr<PileUp>(new PileUp(
        std::move(parts),
        Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
            message.adaptive_step_parameters()),
        Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
            message.fixed_step_parameters()),
        std::move(history),
        psychohistory,
        ephemeris,
        std::move(deletion_callback)));
  }

  if (is_pre_frege) {
    for (auto const& pair : message.actual_part_degrees_of_freedom()) {
      std::uint32_t const part_id = pair.first;
      serialization::Pair const& degrees_of_freedom = pair.second;
      pile_up->actual_part_rigid_motion_.emplace(
          part_id_to_part(part_id),
          RigidMotion<RigidPart, NonRotatingPileUp>::MakeNonRotatingMotion(
              DegreesOfFreedom<NonRotatingPileUp>::ReadFromMessage(
                  degrees_of_freedom)));
    }
    for (auto const& pair : message.apparent_part_degrees_of_freedom()) {
      std::uint32_t const part_id = pair.first;
      serialization::Pair const& degrees_of_freedom = pair.second;
      pile_up->apparent_part_rigid_motion_.emplace(
          part_id_to_part(part_id),
          RigidMotion<RigidPart, ApparentBubble>::MakeNonRotatingMotion(
              DegreesOfFreedom<ApparentBubble>::ReadFromMessage(
                  degrees_of_freedom)));
    }
  } else {
    for (auto const& pair : message.actual_part_rigid_motion()) {
      std::uint32_t const part_id = pair.first;
      serialization::RigidMotion const& rigid_motion = pair.second;
      pile_up->actual_part_rigid_motion_.emplace(
          part_id_to_part(part_id),
          RigidMotion<RigidPart, NonRotatingPileUp>::ReadFromMessage(
              rigid_motion));
    }
    for (auto const& pair : message.apparent_part_rigid_motion()) {
      std::uint32_t const part_id = pair.first;
      serialization::RigidMotion const& rigid_motion = pair.second;
      pile_up->apparent_part_rigid_motion_.emplace(
          part_id_to_part(part_id),
          RigidMotion<RigidPart, ApparentBubble>::ReadFromMessage(
              rigid_motion));
    }
  }
  pile_up->RecomputeFromParts();
  return check_not_null(std::move(pile_up));
}

PileUp::PileUp(
    std::list<not_null<Part*>>&& parts,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters,
    not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> history,
    DiscreteTrajectory<Barycentric>* const psychohistory,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    std::function<void()> deletion_callback)
    : lock_(make_not_null_unique<absl::Mutex>()),
      parts_(std::move(parts)),
      ephemeris_(ephemeris),
      adaptive_step_parameters_(adaptive_step_parameters),
      fixed_step_parameters_(fixed_step_parameters),
      history_(std::move(history)),
      psychohistory_(psychohistory),
      deletion_callback_(std::move(deletion_callback)) {}

void PileUp::DeformPileUpIfNeeded() {
  if (apparent_part_rigid_motion_.empty()) {
    return;
  }
  // A consistency check that |SetPartApparentDegreesOfFreedom| was called for
  // all the parts.
  // TODO(egg): I'd like to log some useful information on check failure, but I
  // need a clean way of getting the debug strings of all parts (rather than
  // giant self-evaluating lambdas).
  CHECK_EQ(parts_.size(), apparent_part_rigid_motion_.size());
  for (not_null<Part*> const part : parts_) {
    CHECK(Contains(apparent_part_rigid_motion_, part));
  }

  // Compute the apparent centre of mass of the parts.
  BarycentreCalculator<DegreesOfFreedom<ApparentBubble>, Mass> calculator;
  for (auto const& pair : apparent_part_rigid_motion_) {
    auto const part = pair.first;
    auto const& apparent_part_rigid_motion = pair.second;
    DegreesOfFreedom<ApparentBubble> const apparent_part_degrees_of_freedom =
        apparent_part_rigid_motion({RigidPart::origin, RigidPart::unmoving});
    calculator.Add(apparent_part_degrees_of_freedom, part->mass());
  }
  auto const apparent_centre_of_mass = calculator.Get();

  // A motion that maps the apparent centre of mass of the parts to the actual
  // centre of mass of the pile-up.
  RigidTransformation<ApparentBubble, NonRotatingPileUp> const
      apparent_bubble_to_pile_up_transformation(
          apparent_centre_of_mass.position(),
          NonRotatingPileUp::origin,
          OrthogonalMap<ApparentBubble, NonRotatingPileUp>::Identity());
  RigidMotion<ApparentBubble, NonRotatingPileUp> const
      apparent_bubble_to_pile_up_motion(
          apparent_bubble_to_pile_up_transformation,
          ApparentBubble::nonrotating,
          apparent_centre_of_mass.velocity());

  // Now update the positions of the parts in the pile-up frame.
  actual_part_rigid_motion_.clear();
  for (auto const& pair : apparent_part_rigid_motion_) {
    auto const part = pair.first;
    auto const& apparent_part_rigid_motion = pair.second;
    actual_part_rigid_motion_.emplace(
        part, apparent_bubble_to_pile_up_motion * apparent_part_rigid_motion);
  }
  apparent_part_rigid_motion_.clear();
}

Status PileUp::AdvanceTime(Instant const& t) {
  CHECK_NOTNULL(psychohistory_);

  Status status;

  Time const Δt = t - psychohistory_->back().time;
  angular_momentum_ += angular_momentum_change_;
  angular_momentum_ += intrinsic_torque_ * Δt;

  auto const history_last = --history_->end();
  if (intrinsic_force_ == Vector<Force, Barycentric>{}) {
    // Remove the fork.
    history_->DeleteFork(psychohistory_);
    if (fixed_instance_ == nullptr) {
      fixed_instance_ = ephemeris_->NewInstance(
          {history_.get()},
          Ephemeris<Barycentric>::NoIntrinsicAccelerations,
          fixed_step_parameters_);
    }
    CHECK_LT(history_->back().time, t);
    status = ephemeris_->FlowWithFixedStep(t, *fixed_instance_);
    psychohistory_ = history_->NewForkAtLast();
    if (history_->back().time < t) {
      // Do not clear the |fixed_instance_| here, we will use it for the next
      // fixed-step integration.
      status.Update(ephemeris_->FlowWithAdaptiveStep(
          psychohistory_,
          Ephemeris<Barycentric>::NoIntrinsicAcceleration,
          t,
          adaptive_step_parameters_,
          Ephemeris<Barycentric>::unlimited_max_ephemeris_steps));
    }
  } else {
    // Destroy the fixed instance, it wouldn't be correct to use it the next
    // time we go through this function.  It will be re-created as needed.
    fixed_instance_ = nullptr;
    // We make the |psychohistory_|, if any, authoritative, i.e. append it to
    // the end of the |history_|. We integrate on top of it, and it gets
    // appended authoritatively to the part tails.
    auto const psychohistory_end = psychohistory_->end();
    auto it = psychohistory_->Fork();
    for (++it; it != psychohistory_end; ++it) {
      history_->Append(it->time, it->degrees_of_freedom);
    }
    history_->DeleteFork(psychohistory_);

    auto const a = intrinsic_force_ / mass_;
    // NOTE(phl): |a| used to be captured by copy below, which is the logical
    // thing to do.  However, since it contains an |R3Element|, it must be
    // aligned on a 16-byte boundary.  Unfortunately, VS2015 gets confused and
    // aligns the function object on an 8-byte boundary, resulting in an
    // addressing fault.  With a reference, VS2015 knows what to do.
    auto const intrinsic_acceleration = [&a](Instant const& t) { return a; };
    status = ephemeris_->FlowWithAdaptiveStep(
        history_.get(),
        intrinsic_acceleration,
        t,
        adaptive_step_parameters_,
        Ephemeris<Barycentric>::unlimited_max_ephemeris_steps);
    psychohistory_ = history_->NewForkAtLast();
  }

  CHECK_NOTNULL(psychohistory_);

  // Append the |history_| authoritatively to the parts' tails and the
  // |psychohistory_| non-authoritatively.
  auto const history_end = history_->end();
  auto const psychohistory_end = psychohistory_->end();
  auto it = history_last;
  for (++it; it != history_end; ++it) {
    AppendToPart<&Part::AppendToHistory>(it);
  }
  it = psychohistory_->Fork();
  for (++it; it != psychohistory_end; ++it) {
    AppendToPart<&Part::AppendToPsychohistory>(it);
  }
  history_->ForgetBefore(psychohistory_->Fork()->time);

  return status;
}

void PileUp::NudgeParts() const {
  auto const actual_centre_of_mass = psychohistory_->back().degrees_of_freedom;

  RigidMotion<Barycentric, NonRotatingPileUp> const barycentric_to_pile_up{
      RigidTransformation<Barycentric, NonRotatingPileUp>{
          actual_centre_of_mass.position(),
          NonRotatingPileUp::origin,
          OrthogonalMap<Barycentric, NonRotatingPileUp>::Identity()},
      Barycentric::nonrotating,
      actual_centre_of_mass.velocity()};
  auto const pile_up_to_barycentric = barycentric_to_pile_up.Inverse();
  for (not_null<Part*> const part : parts_) {
    RigidMotion<RigidPart, Barycentric> const actual_part_rigid_motion =
        pile_up_to_barycentric * FindOrDie(actual_part_rigid_motion_, part);
    part->set_rigid_motion(actual_part_rigid_motion);
  }
}

template<PileUp::AppendToPartTrajectory append_to_part_trajectory>
void PileUp::AppendToPart(DiscreteTrajectory<Barycentric>::Iterator it) const {
  auto const& pile_up_dof = it->degrees_of_freedom;
  RigidMotion<Barycentric, NonRotatingPileUp> const barycentric_to_pile_up(
      RigidTransformation<Barycentric, NonRotatingPileUp>(
          pile_up_dof.position(),
          NonRotatingPileUp::origin,
          OrthogonalMap<Barycentric, NonRotatingPileUp>::Identity()),
      Barycentric::nonrotating,
      pile_up_dof.velocity());
  auto const pile_up_to_barycentric = barycentric_to_pile_up.Inverse();
  for (not_null<Part*> const part : parts_) {
    DegreesOfFreedom<NonRotatingPileUp> const actual_part_degrees_of_freedom =
        FindOrDie(actual_part_rigid_motion_,
                  part)({RigidPart::origin, RigidPart::unmoving});
    (static_cast<Part*>(part)->*append_to_part_trajectory)(
        it->time, pile_up_to_barycentric(actual_part_degrees_of_freedom));
  }
}

PileUpFuture::PileUpFuture(not_null<PileUp const*> const pile_up,
                           std::future<Status> future)
    : pile_up(pile_up), future(std::move(future)) {}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
