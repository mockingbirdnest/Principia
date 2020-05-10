
#include "ksp_plugin/pile_up.hpp"

#include <algorithm>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <utility>

#include "base/flags.hpp"
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
using base::Flags;
using base::make_not_null_unique;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::Commutator;
using geometry::Frame;
using geometry::Identity;
using geometry::NonRotating;
using geometry::Normalize;
using geometry::NormalizeOrZero;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Quaternion;
using geometry::RigidTransformation;
using geometry::Rotation;
using geometry::Velocity;
using geometry::Wedge;
using physics::DegreesOfFreedom;
using physics::RigidMotion;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::AngularMomentum;
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

// NOTE(phl): The default values were tuned on "Fubini oscillator".
constexpr double kp_default = 0.65;
constexpr Inverse<Time> ki_default = 0.03 / Second;
constexpr Time kd_default = 0.0025 * Second;

// Configuration example:
//   principia_flags {
//     kp = 0.65
//     ki = 0.03 / s
//     kd = 0.0025 s
//   }
template<typename Q>
Q GetPIDFlagOr(std::string_view const name, Q const default_value) {
  auto const values = Flags::Values(name);
  if (values.empty()) {
    return default_value;
  } else {
    return ParseQuantity<Q>(*values.cbegin());
  }
}

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
      history_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      deletion_callback_(std::move(deletion_callback)),
      apparent_angular_momentum_controller_(GetPIDFlagOr("kp", kp_default),
                                            GetPIDFlagOr("ki", ki_default),
                                            GetPIDFlagOr("kd", kd_default)) {
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
  MakeEulerSolver(mechanical_system.InertiaTensor(), t);

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
    RigidMotion<RigidPart, Apparent> const& rigid_motion) {
  auto const [_, inserted] =
      apparent_part_rigid_motion_.emplace(part, rigid_motion);
  CHECK(inserted) << "Duplicate part " << part->ShortDebugString() << " at "
                  << rigid_motion;
}

Status PileUp::DeformAndAdvanceTime(Instant const& t) {
  absl::MutexLock l(lock_.get());
  Status status;
  if (psychohistory_->back().time < t) {
    angular_momentum_ += intrinsic_torque_ * (t - psychohistory_->back().time) +
                         angular_momentum_change_;
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
        part_motion.Inverse().angular_velocity_of_to_frame();
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
  history_->WriteToMessage(message->mutable_history(),
                           /*forks=*/{psychohistory_});
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
  bool const is_pre_cesàro = message.history().children().empty();
  bool const is_pre_frege = message.actual_part_degrees_of_freedom_size() > 0 ||
                            message.apparent_part_degrees_of_freedom_size() > 0;
  bool const is_pre_frobenius = message.rigid_pile_up().empty() ||
                                !message.has_angular_momentum();
  std::unique_ptr<PileUp> pile_up;
  if (is_pre_cesàro) {
    if (is_pre_cartan) {
      pile_up = std::unique_ptr<PileUp>(
          new PileUp(std::move(parts),
                     DefaultPsychohistoryParameters(),
                     DefaultHistoryParameters(),
                     DiscreteTrajectory<Barycentric>::ReadFromMessage(
                         message.history(),
                         /*forks=*/{}),
                     /*psychohistory=*/nullptr,
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
              DiscreteTrajectory<Barycentric>::ReadFromMessage(
                  message.history(),
                  /*forks=*/{}),
              /*psychohistory=*/nullptr,
              /*angular_momentum=*/{},
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
    if (is_pre_frobenius) {
      pile_up = std::unique_ptr<PileUp>(
          new PileUp(
              std::move(parts),
              Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
                  message.adaptive_step_parameters()),
              Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
                  message.fixed_step_parameters()),
              std::move(history),
              psychohistory,
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
              std::move(history),
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
              DegreesOfFreedom<Apparent>::ReadFromMessage(
                  degrees_of_freedom)));
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
          RigidMotion<RigidPart, Apparent>::ReadFromMessage(
              rigid_motion));
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
    not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> history,
    DiscreteTrajectory<Barycentric>* const psychohistory,
    Bivector<AngularMomentum, NonRotatingPileUp> const& angular_momentum,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    std::function<void()> deletion_callback)
    : lock_(make_not_null_unique<absl::Mutex>()),
      parts_(std::move(parts)),
      ephemeris_(ephemeris),
      adaptive_step_parameters_(std::move(adaptive_step_parameters)),
      fixed_step_parameters_(std::move(fixed_step_parameters)),
      history_(std::move(history)),
      psychohistory_(psychohistory),
      angular_momentum_(angular_momentum),
      deletion_callback_(std::move(deletion_callback)),
      apparent_angular_momentum_controller_(GetPIDFlagOr("kp", kp_default),
                                            GetPIDFlagOr("ki", ki_default),
                                            GetPIDFlagOr("kd", kd_default)) {}

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
    apparent_angular_momentum_controller_.Clear();
    Bivector<AngularMomentum, PileUpPrincipalAxes> const angular_momentum =
        euler_solver_->AngularMomentumAt(t);
    Rotation<PileUpPrincipalAxes, NonRotatingPileUp> const attitude =
        euler_solver_->AttitudeAt(angular_momentum, t);
    AngularVelocity<NonRotatingPileUp> const angular_velocity_of_rigid_pile_up =
        attitude(euler_solver_->AngularVelocityFor(angular_momentum));

    RigidMotion<PileUpPrincipalAxes, NonRotatingPileUp> const pile_up_motion(
        RigidTransformation<PileUpPrincipalAxes, NonRotatingPileUp>(
            PileUpPrincipalAxes::origin,
            NonRotatingPileUp::origin,
            attitude.Forget<OrthogonalMap>()),
        angular_velocity_of_rigid_pile_up,
        NonRotatingPileUp::unmoving);

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

  auto const angular_momentum_in_apparent_pile_up =
      Identity<NonRotatingPileUp, ApparentPileUp>()(angular_momentum_);

  MechanicalSystem<Apparent, ApparentPileUp> apparent_system;
  for (auto const& [part, apparent_part_rigid_motion] :
       apparent_part_rigid_motion_) {
    apparent_system.AddRigidBody(
        apparent_part_rigid_motion, part->mass(), part->inertia_tensor());
  }
  auto const apparent_centre_of_mass = apparent_system.centre_of_mass();
  auto const apparent_angular_momentum =
      angular_momentum_in_apparent_pile_up -
      apparent_angular_momentum_controller_.ComputeControlVariable(
          apparent_system.AngularMomentum(),
          angular_momentum_in_apparent_pile_up,
          t);
  // Note that the inertia tensor is with respect to the centre of mass, so it
  // is unaffected by the apparent-bubble-to-pile-up correction, which is rigid
  // and involves no change in axes.
  auto const inertia_tensor = apparent_system.InertiaTensor();

  // This is the axis around which we may perform an orientation correction.
  // Identifying the uncorrected (apparent) and corrected coordinates, it is
  // orthogonal to both.
  // Note that the computations are identical in coordinates:
  // |apparent_correction_axis.coordinates() ==
  //  actual_correction_axis.coordinates()|.
  auto const apparent_correction_axis = NormalizeOrZero(Commutator(
      apparent_angular_momentum,
      angular_momentum_in_apparent_pile_up));
  auto const actual_correction_axis = NormalizeOrZero(Commutator(
      Identity<ApparentPileUp, NonRotatingPileUp>()(apparent_angular_momentum),
      angular_momentum_));

  // We apply a rigid rotational correction to the motions of the parts coming
  // from the game (the apparent motions) so as to enforce the conservation of
  // the angular momentum (|angular_momentum_| is authoritative).
  // Mapping L_apparent to L_actual by a rigid motion can be done in many ways.
  // We prefer doing some of that correction by a change of attitude, rather
  // than solely by a change in angular velocity, since, under isotropic
  // conditions, a change in attitude does not alter the physical system (and in
  // particular it does not mess with symplectic integration).  We thus try to
  // map L̂_apparent to L̂_actual, by rotation around the correction axis defined
  // above, which is perpendicular to both, and to impart the appropriate
  // angular velocity correction to correct the norm of L.
  // This amounts to trusting the direction of the angular momentum with respect
  // to the vessel as given to us by the game.
  // However, mapping L̂_apparent to L̂_actual is essentially singular when either
  // is 0.  We remedy to that by only performing the full attitude correction if
  // its angle would be much less than ω Δt, where ω is the geometric mean of
  // the two equivalent angular frequencies |L_actual / I| and |L_apparent / I|.
  // As a result, as either L tends towards 0, so does the largest attitude
  // correction that we allow.

  // The correction is computed via an intermediate frame.
  // In the |EquivalentRigidPileUp| reference frame, a rigid body with the same
  // inertia and angular momentum as the pile up would be immobile.
  // If no attitude correction is performed, the axes of
  // |EquivalentRigidPileUp|, |ApparentPileUp|, and |NonRotatingPileUp| are
  // identical.
  // If an attitude correction is performed, the correction occurs between
  // |ApparentPileUp| and |EquivalentRigidPileUp|.
  // NOTE(egg): while the correction axis is essentially singular at
  // L_apparent = L_actual, the correction itself is only removably singular (as
  // the angle goes to 0).  Since the computation ensured that
  // |apparent_correction_axis.coordinates() ==
  //  actual_correction_axis.coordinates()| exactly, we need not worry about the
  // intermediate essential singularity, and only need to remove it.
  using EquivalentRigidPileUp = Frame<enum class EquivalentRigidPileUpTag>;

  auto const L̂_apparent = NormalizeOrZero(apparent_angular_momentum);
  auto const L̂_actual = NormalizeOrZero(angular_momentum_);
  // NOTE(egg):
  // — The two sides of this disjunction are always equal;
  // — technically this will also be true on the essential singularity.
  bool const on_removable_singularity =
      apparent_correction_axis == Bivector<double, ApparentPileUp>{} ||
      actual_correction_axis == Bivector<double, NonRotatingPileUp>{};
  bool const on_essential_singularity =
      L̂_apparent == Bivector<double, ApparentPileUp>{} ||
      L̂_actual == Bivector<double, NonRotatingPileUp>{};

  bool const trivial_rotations = !correct_orientation ||
                                 on_removable_singularity ||
                                 on_essential_singularity;

  Rotation<ApparentPileUp, EquivalentRigidPileUp>
      apparent_pile_up_equivalent_rotation =
          trivial_rotations
              ? Rotation<ApparentPileUp, EquivalentRigidPileUp>::Identity()
              : Rotation<ApparentPileUp, EquivalentRigidPileUp>(
                    apparent_correction_axis,
                    L̂_apparent,
                    Commutator(apparent_correction_axis, L̂_apparent));
  Rotation<NonRotatingPileUp, EquivalentRigidPileUp>
      actual_pile_up_equivalent_rotation =
          trivial_rotations
              ? Rotation<NonRotatingPileUp, EquivalentRigidPileUp>::Identity()
              : Rotation<NonRotatingPileUp, EquivalentRigidPileUp>(
                    actual_correction_axis,
                    L̂_actual,
                    Commutator(actual_correction_axis, L̂_actual));
  Rotation<ApparentPileUp, NonRotatingPileUp> const
      tentative_attitude_correction =
          actual_pile_up_equivalent_rotation.Inverse() *
          apparent_pile_up_equivalent_rotation;

  // The angular velocity of a rigid body with the inertia and angular momentum
  // of the apparent parts.
  auto const apparent_equivalent_angular_velocity =
      apparent_angular_momentum / inertia_tensor;
  // The angular velocity of a rigid body with the inertia of the apparent
  // parts, and the angular momentum of the pile up.
  auto actual_equivalent_angular_velocity =
      angular_momentum_ / tentative_attitude_correction(inertia_tensor);

  // So far we have only dealt with the essential singularity by removing it.
  // We now need to deal with its neighbourhood, by ensuring that the tentative
  // correction is not too big compared to the angular velocities involved.
  Quaternion const q = tentative_attitude_correction.quaternion();
  // α is the angle of the tentative attitude correction.
  Angle const α =
      2 * quantities::ArcTan(q.imaginary_part().Norm(), q.real_part());
  // The geometric mean of the norms of the equivalent angular velocities; this
  // is zero when either one is zero, and thus it is zero when either angular
  // momentum is zero.
  AngularFrequency const ω =
      thresholding ? std::min(apparent_equivalent_angular_velocity.Norm(),
                              actual_equivalent_angular_velocity.Norm())
                   : Sqrt(apparent_equivalent_angular_velocity.Norm() *
                          actual_equivalent_angular_velocity.Norm());
  Time const Δt = t - psychohistory_->back().time;
  // This correction angle is bounded by ω Δt (correct at most as fast as the
  // ship is turning) and by α (correct at most what is needed).
  // This in particular eliminates the neighbourhood of the singularity where L
  // is small (and thus ω is small), by using α itself as a scale-free
  // definition of “small”.
  Angle const β = thresholding ? (α > ω * Δt ? Angle{} : α)
                               : Tanh(ω * Δt / α * Radian) * α;
  auto attitude_correction =
      Rotation<ApparentPileUp, EquivalentRigidPileUp>::Identity();
  if (!trivial_rotations) {
    // An active rotation of angle β around apparent_correction_axis maps
    // L̂_apparent to its corrected direction.
    // TODO(egg): The definitions might be more readable if
    // |EquivalentRigidPileUp| were defined as in Fubini, with its y axis being
    // L̂_apparent in |ApparentPileUp| and
    //   Identity<ApparentPileUp, NonRotatingPileUp>()(
    //       Rotation<ApparentPileUp, ApparentPileUp>(β)(L̂_apparent));
    // in |NonRotatingPileUp|.  On the other hand, this would make the
    // computations even more convoluted.
    attitude_correction =
        Rotation<ApparentPileUp, EquivalentRigidPileUp>::Identity() *
        Rotation<ApparentPileUp, ApparentPileUp>(β, apparent_correction_axis);
  }
  actual_equivalent_angular_velocity =
      angular_momentum_ / Identity<EquivalentRigidPileUp, NonRotatingPileUp>()(
                              attitude_correction(inertia_tensor));

  bool const correcting_orientation =
      apparent_pile_up_equivalent_rotation.quaternion() != Quaternion(1) ||
      actual_pile_up_equivalent_rotation.quaternion() != Quaternion(1);

  RigidMotion<ApparentPileUp, EquivalentRigidPileUp> const
      apparent_pile_up_equivalent_motion(
          RigidTransformation<ApparentPileUp, EquivalentRigidPileUp>(
              ApparentPileUp::origin,
              EquivalentRigidPileUp::origin,
              attitude_correction.Forget<OrthogonalMap>()),
          correct_angular_velocity ? apparent_equivalent_angular_velocity
                                   : ApparentPileUp::nonrotating,
          ApparentPileUp::unmoving);
  RigidMotion<NonRotatingPileUp, EquivalentRigidPileUp> const
      actual_pile_up_equivalent_motion(
          RigidTransformation<NonRotatingPileUp, EquivalentRigidPileUp>(
              NonRotatingPileUp::origin,
              EquivalentRigidPileUp::origin,
              OrthogonalMap<NonRotatingPileUp,
                            EquivalentRigidPileUp>::Identity()),
          correct_angular_velocity ? actual_equivalent_angular_velocity
                                   : NonRotatingPileUp::nonrotating,
          NonRotatingPileUp::unmoving);
  RigidMotion<Apparent, NonRotatingPileUp> const
      apparent_to_pile_up_motion =
          actual_pile_up_equivalent_motion.Inverse() *
          apparent_pile_up_equivalent_motion *
          apparent_system.LinearMotion().Inverse();

  // Now update the motions of the parts in the pile-up frame.
  actual_part_rigid_motion_.clear();
  for (auto const& pair : apparent_part_rigid_motion_) {
    auto const part = pair.first;
    auto const& apparent_part_rigid_motion = pair.second;
    actual_part_rigid_motion_.emplace(
        part, apparent_to_pile_up_motion * apparent_part_rigid_motion);
  }
  apparent_part_rigid_motion_.clear();

  std::stringstream s;
  s << "Apparent: " << apparent_angular_momentum << "\n"
    << "norm: " << apparent_angular_momentum.Norm() << "\n"
    << "Actual: " << angular_momentum_ << "\n"
    << "norm: " << angular_momentum_.Norm() << "\n"
    << "|Lap-Lac|: "
    << (angular_momentum_ - Identity<ApparentPileUp, NonRotatingPileUp>()(
                                apparent_angular_momentum))
           .Norm() << "\n"
    << "|Lap|-|Lac|: "
    << angular_momentum_.Norm() - apparent_angular_momentum.Norm() << "\n"
    << u8"∡Lap, Lac: "
    << geometry::AngleBetween(angular_momentum_,
                              Identity<ApparentPileUp, NonRotatingPileUp>()(
                                  apparent_angular_momentum)) /
           quantities::si::Degree
    << u8"°\n"
    << u8"α: " << α / quantities::si::Degree << u8"°\n"
    << u8"ωΔt: " << ω * Δt / quantities::si::Degree << u8"°\n"
    << u8"β: " << β / quantities::si::Degree << u8"°\n"
    << u8"ω: " << ω / (2 * π * Radian / quantities::si::Minute) << " rpm\n"
    << u8"|ωap|: "
    << apparent_equivalent_angular_velocity.Norm() /
           (2 * π * Radian / quantities::si::Minute)
    << " rpm\n"
    << u8"|ωac|: "
    << actual_equivalent_angular_velocity.Norm() /
           (2 * π * Radian / quantities::si::Minute)
    << " rpm\n"
    << u8"Δt: " << Δt / quantities::si::Milli(quantities::si::Second)
    << u8"ms\n";
  trace = s.str();

  MakeEulerSolver(Identity<ApparentPileUp, NonRotatingPileUp>()(inertia_tensor),
                  t);
}

Status PileUp::AdvanceTime(Instant const& t) {
  CHECK_NOTNULL(psychohistory_);

  Status status;
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
      status.Update(
          ephemeris_->FlowWithAdaptiveStep(
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

  // Append the |history_| to the parts' history and the |psychohistory_| to the
  // parts' psychohistory.  Drop the history of the pile-up, we won't need it
  // anymore.
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
        FindOrDie(actual_part_rigid_motion_, part)({RigidPart::origin,
                                                    RigidPart::unmoving});
    (static_cast<Part*>(part)->*append_to_part_trajectory)(
        it->time,
        pile_up_to_barycentric(actual_part_degrees_of_freedom));
  }
}

PileUpFuture::PileUpFuture(not_null<PileUp const*> const pile_up,
                           std::future<Status> future)
    : pile_up(pile_up),
      future(std::move(future)) {}

bool PileUp::correct_orientation = true;
bool PileUp::correct_angular_velocity = true;
bool PileUp::thresholding = false;

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
