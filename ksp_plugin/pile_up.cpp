
#include "ksp_plugin/pile_up.hpp"

#include <algorithm>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <utility>

#include "base/map_util.hpp"
#include "geometry/identity.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/part.hpp"
#include "quantities/parser.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

using base::check_not_null;
using base::FindOrDie;
using base::make_not_null_unique;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::Commutator;
using geometry::DeduceSignPreservingOrientation;
using geometry::DeduceSignReversingOrientation;
using geometry::EvenPermutation;
using geometry::Frame;
using geometry::Identity;
using geometry::NonRotating;
using geometry::Normalize;
using geometry::NormalizeOrZero;
using geometry::OddPermutation;
using geometry::OrthogonalMap;
using geometry::Permutation;
using geometry::Position;
using geometry::Quaternion;
using geometry::RigidTransformation;
using geometry::Rotation;
using geometry::Sign;
using geometry::Signature;
using geometry::Velocity;
using geometry::Wedge;
using physics::DegreesOfFreedom;
using physics::RigidMotion;
using quantities::Abs;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::AngularMomentum;
using quantities::Infinity;
using quantities::Inverse;
using quantities::ParseQuantity;
using quantities::Sqrt;
using quantities::Tanh;
using quantities::Time;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

PileUp::PileUp(
    std::list<not_null<Part*>>&& parts,
    Instant const& t,
    Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters,
    Ephemeris<Barycentric>::FixedStepParameters fixed_step_parameters,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    std::function<void()> deletion_callback)
    : lock_(make_not_null_unique<absl::Mutex>()),
      parts_(std::move(parts)),
      ephemeris_(ephemeris),
      adaptive_step_parameters_(std::move(adaptive_step_parameters)),
      fixed_step_parameters_(std::move(fixed_step_parameters)),
      history_(trajectory_.segments().begin()),
      deletion_callback_(std::move(deletion_callback)) {
  LOG(INFO) << "Constructing pile up at " << this;
  MechanicalSystem<Barycentric, NonRotatingPileUp> mechanical_system;
  for (not_null<Part*> const part : parts_) {
    mechanical_system.AddRigidBody(
        part->rigid_motion(), part->mass(), part->inertia_tensor());
  }
  auto const barycentre = mechanical_system.centre_of_mass();
  trajectory_.Append(t, barycentre);

  angular_momentum_ = mechanical_system.AngularMomentum();

  RigidMotion<Barycentric, NonRotatingPileUp> const barycentric_to_pile_up =
      mechanical_system.LinearMotion().Inverse();
  for (not_null<Part*> const part : parts_) {
    actual_part_rigid_motion_.emplace(
        part, barycentric_to_pile_up * part->rigid_motion());
  }
  MakeEulerSolver(mechanical_system.InertiaTensor(), t);

  psychohistory_ = trajectory_.NewSegment();

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
    RigidMotion<RigidPart, Apparent> const& rigid_motion) {
  auto const [_, inserted] =
      apparent_part_rigid_motion_.emplace(part, rigid_motion);
  CHECK(inserted) << "Duplicate part " << part->ShortDebugString() << " at "
                  << rigid_motion;
}

absl::Status PileUp::DeformAndAdvanceTime(Instant const& t) {
  absl::MutexLock l(lock_.get());
  absl::Status status;
  if (psychohistory_->back().time < t) {
    DeformPileUpIfNeeded(t);
    status = AdvanceTime(t);
    NudgeParts();
  }
  return status;
}

void PileUp::RecomputeFromParts() {
  absl::MutexLock l(lock_.get());
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
        part_motion.angular_velocity_of<RigidPart>();
    InertiaTensor<NonRotatingPileUp> part_inertia_tensor =
        part_motion.orthogonal_map()(part->inertia_tensor());
    if (part->is_solid_rocket_motor()) {
      // KSP makes the inertia tensor vary proportionally to the mass; this
      // corresponds to the body uniformly changing density.
      angular_momentum_change_ +=
          Wedge(part_dof.position() - NonRotatingPileUp::origin,
                part->mass_change() * part_dof.velocity()) *
              Radian +
          part->mass_change() / part->mass() *
              (part_inertia_tensor * part_angular_velocity);
    }
  }
}

void PileUp::WriteToMessage(not_null<serialization::PileUp*> message) const {
  for (not_null<Part*> const part : parts_) {
    message->add_part_id(part->part_id());
  }
  trajectory_.WriteToMessage(message->mutable_history(),
                             /*forks=*/{history_, psychohistory_},
                             /*exact=*/{});
  for (auto const& [part, rigid_motion] : actual_part_rigid_motion_) {
    rigid_motion.WriteToMessage(&(
        (*message->mutable_actual_part_rigid_motion())[part->part_id()]));
  }
  for (auto const& [part, rigid_motion] : apparent_part_rigid_motion_) {
    rigid_motion.WriteToMessage(&(
        (*message->mutable_apparent_part_rigid_motion())[part->part_id()]));
  }
  for (auto const& [part, rigid_transformation] : rigid_pile_up_) {
    rigid_transformation.WriteToMessage(&(
        (*message->mutable_rigid_pile_up())[part->part_id()]));
  }
  if (euler_solver_.has_value()) {
    euler_solver_->WriteToMessage(message->mutable_euler_solver());
  }
  angular_momentum_.WriteToMessage(message->mutable_angular_momentum());
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
  bool const is_pre_cesàro = message.history().children().empty() &&
                             message.history().segment_size() == 0;
  bool const is_pre_frege = message.actual_part_degrees_of_freedom_size() > 0 ||
                            message.apparent_part_degrees_of_freedom_size() > 0;
  bool const is_pre_frobenius = message.rigid_pile_up().empty() ||
                                !message.has_angular_momentum();
  bool const is_pre_ζήνων = message.history().segment_size() == 0;
  LOG_IF(WARNING, is_pre_frobenius)
      << "Reading pre-"
      << (is_pre_cartan   ? "Cartan"
          : is_pre_cesàro ? u8"Cesàro"
          : is_pre_frege  ? "Frege"
          : is_pre_ζήνων  ? "Frobenius"
                          : u8"pre-Ζήνων") << " PileUp";

  std::unique_ptr<PileUp> pile_up;
  if (is_pre_cesàro) {
    if (is_pre_cartan) {
      pile_up = std::unique_ptr<PileUp>(
          new PileUp(std::move(parts),
                     DefaultPsychohistoryParameters(),
                     DefaultHistoryParameters(),
                     DiscreteTraject0ry<Barycentric>::ReadFromMessage(
                         message.history(),
                         /*forks=*/{}),
                     /*history=*/std::nullopt,
                     /*psychohistory=*/std::nullopt,
                     /*angular_momentum=*/{},
                     ephemeris,
                     std::move(deletion_callback)));
    } else {
      pile_up = std::unique_ptr<PileUp>(
          new PileUp(
              std::move(parts),
              Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
                  message.adaptive_step_parameters()),
              Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
                  message.fixed_step_parameters()),
              DiscreteTraject0ry<Barycentric>::ReadFromMessage(
                  message.history(),
                  /*forks=*/{}),
              /*history=*/std::nullopt,
              /*psychohistory=*/std::nullopt,
              /*angular_momentum=*/{},
              ephemeris,
              std::move(deletion_callback)));
    }
    // Fork a psychohistory for compatibility if there is a non-authoritative
    // point.
    if (pile_up->history_->size() == 2) {
      DiscreteTraject0ry<Barycentric> psychohistory;
      auto const history_second_point = std::next(pile_up->history_->begin());
      for (auto it = history_second_point;
           it != pile_up->history_->end();
           ++it) {
        psychohistory.Append(it->time, it->degrees_of_freedom);
      }
      pile_up->trajectory_.ForgetAfter(history_second_point);
      pile_up->psychohistory_ =
          pile_up->trajectory_.AttachSegments(std::move(psychohistory));
    } else {
      pile_up->psychohistory_ = pile_up->trajectory_.NewSegment();
    }
  } else {
    if (is_pre_frobenius) {
      DiscreteTrajectorySegmentIterator<Barycentric> psychohistory;
      auto trajectory = DiscreteTraject0ry<Barycentric>::ReadFromMessage(
          message.history(),
          /*forks=*/{&psychohistory});
      pile_up = std::unique_ptr<PileUp>(
          new PileUp(
              std::move(parts),
              Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
                  message.adaptive_step_parameters()),
              Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
                  message.fixed_step_parameters()),
              std::move(trajectory),
              /*history=*/std::nullopt,
              psychohistory,
              /*angular_momentum=*/{},
              ephemeris,
              std::move(deletion_callback)));
    } else if (is_pre_ζήνων) {
      DiscreteTrajectorySegmentIterator<Barycentric> psychohistory;
      auto trajectory = DiscreteTraject0ry<Barycentric>::ReadFromMessage(
          message.history(),
          /*forks=*/{&psychohistory});
      pile_up = std::unique_ptr<PileUp>(
          new PileUp(
              std::move(parts),
              Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
                  message.adaptive_step_parameters()),
              Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
                  message.fixed_step_parameters()),
              std::move(trajectory),
              /*history=*/std::nullopt,
              psychohistory,
              Bivector<AngularMomentum, NonRotatingPileUp>::ReadFromMessage(
                  message.angular_momentum()),
              ephemeris,
              std::move(deletion_callback)));
    } else {
      DiscreteTrajectorySegmentIterator<Barycentric> history;
      DiscreteTrajectorySegmentIterator<Barycentric> psychohistory;
      auto trajectory = DiscreteTraject0ry<Barycentric>::ReadFromMessage(
          message.history(),
          /*forks=*/{&history, &psychohistory});
      pile_up = std::unique_ptr<PileUp>(
          new PileUp(
              std::move(parts),
              Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
                  message.adaptive_step_parameters()),
              Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
                  message.fixed_step_parameters()),
              std::move(trajectory),
              history,
              psychohistory,
              Bivector<AngularMomentum, NonRotatingPileUp>::ReadFromMessage(
                  message.angular_momentum()),
              ephemeris,
              std::move(deletion_callback)));
    }
  }

  if (is_pre_frege) {
    for (auto const& [part_id, degrees_of_freedom] :
         message.actual_part_degrees_of_freedom()) {
      pile_up->actual_part_rigid_motion_.emplace(
          part_id_to_part(part_id),
          RigidMotion<RigidPart, NonRotatingPileUp>::MakeNonRotatingMotion(
              DegreesOfFreedom<NonRotatingPileUp>::ReadFromMessage(
                  degrees_of_freedom)));
    }
    for (auto const& [part_id, degrees_of_freedom] :
         message.apparent_part_degrees_of_freedom()) {
      pile_up->apparent_part_rigid_motion_.emplace(
          part_id_to_part(part_id),
          RigidMotion<RigidPart, Apparent>::MakeNonRotatingMotion(
              DegreesOfFreedom<Apparent>::ReadFromMessage(degrees_of_freedom)));
    }
  } else {
    for (auto const& [part_id, rigid_motion] :
         message.actual_part_rigid_motion()) {
      pile_up->actual_part_rigid_motion_.emplace(
          part_id_to_part(part_id),
          RigidMotion<RigidPart, NonRotatingPileUp>::ReadFromMessage(
              rigid_motion));
    }
    for (auto const& [part_id, rigid_motion] :
         message.apparent_part_rigid_motion()) {
      pile_up->apparent_part_rigid_motion_.emplace(
          part_id_to_part(part_id),
          RigidMotion<RigidPart, Apparent>::ReadFromMessage(rigid_motion));
    }
  }

  if (is_pre_frobenius) {
    MechanicalSystem<Barycentric, NonRotatingPileUp> mechanical_system;
    for (not_null<Part*> const part : pile_up->parts_) {
      mechanical_system.AddRigidBody(
          part->rigid_motion(), part->mass(), part->inertia_tensor());
    }
    pile_up->MakeEulerSolver(mechanical_system.InertiaTensor(),
                             pile_up->psychohistory_->back().time);
  } else {
    for (auto const& [part_id, rigid_transformation] :
         message.rigid_pile_up()) {
      pile_up->rigid_pile_up_.emplace(
          part_id_to_part(part_id),
          RigidTransformation<RigidPart, PileUpPrincipalAxes>::ReadFromMessage(
              rigid_transformation));
    }
    if (message.has_euler_solver()) {
      pile_up->euler_solver_.emplace(
          EulerSolver<NonRotatingPileUp, PileUpPrincipalAxes>::ReadFromMessage(
              message.euler_solver()));
    }
  }
  pile_up->RecomputeFromParts();

  return check_not_null(std::move(pile_up));
}

PileUp::PileUp(
    std::list<not_null<Part*>>&& parts,
    Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters,
    Ephemeris<Barycentric>::FixedStepParameters fixed_step_parameters,
    DiscreteTraject0ry<Barycentric> trajectory,
    std::optional<DiscreteTrajectorySegmentIterator<Barycentric>> history,
    std::optional<DiscreteTrajectorySegmentIterator<Barycentric>> psychohistory,
    Bivector<AngularMomentum, NonRotatingPileUp> const& angular_momentum,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    std::function<void()> deletion_callback)
    : lock_(make_not_null_unique<absl::Mutex>()),
      parts_(std::move(parts)),
      ephemeris_(ephemeris),
      adaptive_step_parameters_(std::move(adaptive_step_parameters)),
      fixed_step_parameters_(std::move(fixed_step_parameters)),
      trajectory_(std::move(trajectory)),
      angular_momentum_(angular_momentum),
      deletion_callback_(std::move(deletion_callback)) {
  if (history.has_value()) {
    history_ = history.value();
  } else {
    history_ = trajectory_.segments().begin();
  }
  if (psychohistory.has_value()) {
    psychohistory_ = psychohistory.value();
  }
}

void PileUp::MakeEulerSolver(
    InertiaTensor<NonRotatingPileUp> const& inertia_tensor,
    Instant const& t) {
  auto const eigensystem = inertia_tensor.Diagonalize<PileUpPrincipalAxes>();
  euler_solver_.emplace(eigensystem.form.coordinates().Diagonal(),
                        angular_momentum_,
                        eigensystem.rotation,
                        t);
  RigidTransformation<NonRotatingPileUp, PileUpPrincipalAxes> const
      to_pile_up_principal_axes(
          NonRotatingPileUp::origin,
          PileUpPrincipalAxes::origin,
          eigensystem.rotation.Inverse().Forget<OrthogonalMap>());
  rigid_pile_up_.clear();
  for (auto const& [part, actual_rigid_motion] : actual_part_rigid_motion_) {
    rigid_pile_up_.emplace(
        part,
        to_pile_up_principal_axes * actual_rigid_motion.rigid_transformation());
  }
}

void PileUp::DeformPileUpIfNeeded(Instant const& t) {
  if (apparent_part_rigid_motion_.empty()) {
    RigidMotion<PileUpPrincipalAxes, NonRotatingPileUp> const pile_up_motion =
        euler_solver_->MotionAt(
            t, {NonRotatingPileUp::origin, NonRotatingPileUp::unmoving});

    for (auto& [part, actual_rigid_motion] : actual_part_rigid_motion_) {
      actual_rigid_motion =
          pile_up_motion * RigidMotion<RigidPart, PileUpPrincipalAxes>(
                               rigid_pile_up_.at(part),
                               PileUpPrincipalAxes::nonrotating,
                               PileUpPrincipalAxes::unmoving);
    }
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

  Instant const& t0 = psychohistory_->back().time;
  Time const Δt = t - t0;

  MechanicalSystem<Apparent, ApparentPileUp> apparent_system;
  for (auto const& [part, apparent_part_rigid_motion] :
       apparent_part_rigid_motion_) {
    apparent_system.AddRigidBody(
        apparent_part_rigid_motion, part->mass(), part->inertia_tensor());
  }
  auto const apparent_angular_momentum = apparent_system.AngularMomentum();
  auto const apparent_inertia_tensor = apparent_system.InertiaTensor();
  auto const apparent_inertia_eigensystem =
      apparent_inertia_tensor.Diagonalize<PileUpPrincipalAxes>();

  Rotation<PileUpPrincipalAxes, ApparentPileUp> const apparent_attitude =
      apparent_inertia_eigensystem.rotation;

  // The motion of a hypothetical rigid body with the same moment of inertia and
  // angular momentum as the apparent parts.
  RigidMotion<PileUpPrincipalAxes, ApparentPileUp> const
      apparent_pile_up_motion(
          RigidTransformation<PileUpPrincipalAxes, ApparentPileUp>(
              PileUpPrincipalAxes::origin,
              ApparentPileUp::origin,
              apparent_attitude.Forget<OrthogonalMap>()),
          apparent_angular_momentum / apparent_inertia_tensor,
          ApparentPileUp::unmoving);

  // In a non-rigid body, the principal axes are not stable, and cannot be used
  // to determine attitude.  We treat this as a flexible body, and use a part to
  // propagate the attitude: its attitude at |t0| is used, together with its
  // orientation with respect to the new principal axes (from apparent
  // coordinates at |t|), to define the attitude at |t0| of the new principal
  // axes.
  // We choose the part which rotates the least with respect to the
  // (instantaneous) principal axes to define the attitude.
  AngularFrequency reference_part_proper_ω = Infinity<AngularFrequency>;
  Part* reference_part = nullptr;
  std::optional<OrthogonalMap<RigidPart, PileUpPrincipalAxes>>
      reference_part_orientation_in_principal_axes;
  for (auto const& [part, apparent_part_motion] : apparent_part_rigid_motion_) {
    RigidMotion<RigidPart, PileUpPrincipalAxes> const
        part_motion_in_principal_axes =
            apparent_pile_up_motion.Inverse() *
            apparent_system.LinearMotion().Inverse() * apparent_part_motion;
    auto const part_proper_ω =
        part_motion_in_principal_axes.angular_velocity_of<PileUpPrincipalAxes>()
            .Norm();
    if (part_proper_ω < reference_part_proper_ω) {
      reference_part_proper_ω = part_proper_ω;
      reference_part = part;
      reference_part_orientation_in_principal_axes =
          part_motion_in_principal_axes.orthogonal_map();
    }
  }

  OrthogonalMap<RigidPart, NonRotatingPileUp> const
      reference_part_initial_attitude =
          actual_part_rigid_motion_.at(reference_part).orthogonal_map();

  // |reference_part_initial_attitude|, an |actual_part_rigid_motion_|, is
  // computed from the output of the previous |EulerSolver|.  |initial_attitude|
  // will be used to construct the new one.  In order to prevent roundoff
  // accumulation from eventually producing noticeably non-unit quaternions, we
  // normalize |initial_attitude|.
  Rotation<PileUpPrincipalAxes, NonRotatingPileUp> initial_attitude =
      (reference_part_initial_attitude *
       reference_part_orientation_in_principal_axes->Inverse())
          .AsRotation();
  initial_attitude = Rotation<PileUpPrincipalAxes, NonRotatingPileUp>(
      Normalize(initial_attitude.quaternion()));

  // We take into account the changes to |angular_momentum_| and to the moments
  // of inertia for the step from t0 to t before propagating the attitude from
  // t0 to t. This forms a splitting with the game, with the game changing
  // angular momentum and moment of inertia according to various physical
  // effects (engines, aerodynamics, internal dynamics, etc.) that depend, among
  // other things, on the attitude and angular velocity, and the Euler solver
  // changing attitude and angular velocity according to Euler’s equations.
  angular_momentum_ += intrinsic_torque_ * Δt + angular_momentum_change_;
  euler_solver_.emplace(
      apparent_inertia_eigensystem.form.coordinates().Diagonal(),
      angular_momentum_,
      initial_attitude,
      t0);

  // This is where we compute our half of the splitting.
  RigidMotion<PileUpPrincipalAxes, NonRotatingPileUp> const
      actual_pile_up_motion = euler_solver_->MotionAt(
          t, {NonRotatingPileUp::origin, NonRotatingPileUp::unmoving});

  RigidMotion<ApparentPileUp, NonRotatingPileUp> const rotational_correction =
      actual_pile_up_motion * apparent_pile_up_motion.Inverse();
  RigidMotion<Apparent, NonRotatingPileUp> const correction =
      rotational_correction * apparent_system.LinearMotion().Inverse();

  // Now update the motions of the parts in the pile-up frame, and keep their
  // orientations with respect to the principal axes in case we warp.
  actual_part_rigid_motion_.clear();
  rigid_pile_up_.clear();
  for (auto const& [part, apparent_part_rigid_motion] :
       apparent_part_rigid_motion_) {
    RigidMotion<RigidPart, NonRotatingPileUp> const actual_rigid_motion =
        correction * apparent_part_rigid_motion;

    actual_part_rigid_motion_.emplace(part, actual_rigid_motion);
    rigid_pile_up_.emplace(
        part,
        actual_pile_up_motion.rigid_transformation().Inverse() *
            actual_rigid_motion.rigid_transformation());
  }
  apparent_part_rigid_motion_.clear();
}

absl::Status PileUp::AdvanceTime(Instant const& t) {
  absl::Status status;
  Instant const history_last = history_->back().time;
  if (intrinsic_force_ == Vector<Force, Barycentric>{}) {
    // Remove the fork.
    trajectory_.DeleteSegments(psychohistory_);
    if (fixed_instance_ == nullptr) {
      fixed_instance_ = ephemeris_->NewInstance(
          {&trajectory_},
          Ephemeris<Barycentric>::NoIntrinsicAccelerations,
          fixed_step_parameters_);
    }
    CHECK_LT(history_->back().time, t);
    status = ephemeris_->FlowWithFixedStep(t, *fixed_instance_);
    psychohistory_ = trajectory_.NewSegment();
    if (history_->back().time < t) {
      // Do not clear the |fixed_instance_| here, we will use it for the next
      // fixed-step integration.
      status.Update(
          ephemeris_->FlowWithAdaptiveStep(
              &trajectory_,
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
    auto const psychohistory_trajectory =
        trajectory_.DetachSegments(psychohistory_);
    for (auto const& [time, degrees_of_freedom] : psychohistory_trajectory) {
      trajectory_.Append(time, degrees_of_freedom);
    }

    auto const intrinsic_acceleration =
        [a = intrinsic_force_ / mass_](Instant const& t) { return a; };
    status = ephemeris_->FlowWithAdaptiveStep(
                 &trajectory_,
                 intrinsic_acceleration,
                 t,
                 adaptive_step_parameters_,
                 Ephemeris<Barycentric>::unlimited_max_ephemeris_steps);
    psychohistory_ = trajectory_.NewSegment();
  }

  // Append the |history_| to the parts' history and the |psychohistory_| to the
  // parts' psychohistory.  Drop the history of the pile-up, we won't need it
  // anymore.
  auto const history_end = history_->end();
  auto const psychohistory_end = psychohistory_->end();
  for (auto it = trajectory_.upper_bound(history_last);
       it != history_end;
       ++it) {
    AppendToPart<&Part::AppendToHistory>(it);
  }
  for (auto it = history_end; it != psychohistory_end; ++it) {
    AppendToPart<&Part::AppendToPsychohistory>(it);
  }
  trajectory_.ForgetBefore(psychohistory_->front().time);

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
void PileUp::AppendToPart(DiscreteTraject0ry<Barycentric>::iterator it) const {
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
        FindOrDie(actual_part_rigid_motion_, part)({RigidPart::origin,
                                                    RigidPart::unmoving});
    (static_cast<Part*>(part)->*append_to_part_trajectory)(
        it->time,
        pile_up_to_barycentric(actual_part_degrees_of_freedom));
  }
}

PileUpFuture::PileUpFuture(not_null<PileUp const*> const pile_up,
                           std::future<absl::Status> future)
    : pile_up(pile_up),
      future(std::move(future)) {}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
