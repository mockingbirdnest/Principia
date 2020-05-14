
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
using geometry::DeduceSignPreservingOrientation;
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
                                            GetPIDFlagOr("kd", kd_default)),
      logger_(SOLUTION_DIR / "mathematica" / "pile_up", true) {
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
    AdvanceEulerSolver(t);
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
                                            GetPIDFlagOr("kd", kd_default)),
      logger_(SOLUTION_DIR / "mathematica" / "pile_up", true) {}

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

void PileUp::AdvanceEulerSolver(Instant t) {
  auto const t0 = psychohistory_->back().time;
  auto const Δt = t - t0;
  // TODO(egg): We integrate the |angular_momentum_change_| resulting from
  // changes in mass, but we do not take into account the change in inertia
  // tensor here.
  auto const effective_torque =
      intrinsic_torque_ + angular_momentum_change_ / Δt;
  if (effective_torque == Bivector<Torque, NonRotatingPileUp>{}) {
    return;
  }
  Rotation<PileUpPrincipalAxes, NonRotatingPileUp> const initial_attitude =
      euler_solver_->AttitudeAt(euler_solver_->AngularMomentumAt(t0), t0);
  if (instant_initial_forces) {
    angular_momentum_ += intrinsic_torque_ * Δt + angular_momentum_change_;
    euler_solver_.emplace(euler_solver_->moments_of_inertia(),
                          angular_momentum_,
                          initial_attitude,
                          t0);
  } else {
    Bivector<Torque, PileUpPrincipalAxes> effective_torque_in_principal_axes =
        initial_attitude.Inverse()(effective_torque);
    // A hundred steps of NewtonDelambreStørmerVerletLeapfrog at 50 Hz, more in
    // physics warp.  This is overkill, but it makes sure we are integrating
    // correctly for all practical purposes.
    // TODO(egg): integrators for |DecomposableFirstOrderDifferentialEquation|.
    int const n = std::round(Δt / (200 * quantities::si::Micro(Second)));
    Time const h = Δt / n;
    for (auto [i, ti] = std::pair{0, t0}; i < n; ++i, ti = t0 + i * h) {
      Rotation<PileUpPrincipalAxes, NonRotatingPileUp> attitude_at_first_stage =
          euler_solver_->AttitudeAt(
              euler_solver_->AngularMomentumAt(ti + h / 2), ti + h / 2);
      // TODO(phl): It seems that the solver doesn’t always provide normalized
      // quaternions; this adds up over the steps and causes the rocket to blow
      // up.
      attitude_at_first_stage =
          Rotation<PileUpPrincipalAxes, NonRotatingPileUp>(
              attitude_at_first_stage.quaternion() /
              attitude_at_first_stage.quaternion().Norm());
      if (inertially_fixed_forces) {
        for (not_null<Part*> const part : parts_) {
          angular_momentum_ +=
              h * (Wedge(attitude_at_first_stage(
                             rigid_pile_up_.at(part)(RigidPart::origin) -
                             PileUpPrincipalAxes::origin),
                         Identity<Barycentric, NonRotatingPileUp>()(
                             part->intrinsic_force())) *
                       Radian +
                   Identity<Barycentric, NonRotatingPileUp>()(
                       part->intrinsic_torque()));
        }
      } else {
        angular_momentum_ += body_fixed_forces
                                 ? attitude_at_first_stage(
                                       effective_torque_in_principal_axes * h)
                                 : effective_torque * h;
      }
      euler_solver_.emplace(euler_solver_->moments_of_inertia(),
                            angular_momentum_,
                            attitude_at_first_stage,
                            ti + h / 2);
    }
  }
}

void PileUp::DeformPileUpIfNeeded(Instant const& t) {
  // This is the motion that the pile up would have if it were rigid and had
  // been subjected to constant intrinsic torque and angular momentum change.
  RigidMotion<PileUpPrincipalAxes,
              NonRotatingPileUp> const pile_up_actual_motion = [this, t] {
    Bivector<AngularMomentum, PileUpPrincipalAxes> const angular_momentum =
        euler_solver_->AngularMomentumAt(t);
    Rotation<PileUpPrincipalAxes, NonRotatingPileUp> const attitude =
        euler_solver_->AttitudeAt(angular_momentum, t);
    AngularVelocity<NonRotatingPileUp> const angular_velocity_of_rigid_pile_up =
        attitude(euler_solver_->AngularVelocityFor(angular_momentum));

    return RigidMotion<PileUpPrincipalAxes, NonRotatingPileUp>(
        RigidTransformation<PileUpPrincipalAxes, NonRotatingPileUp>(
            PileUpPrincipalAxes::origin,
            NonRotatingPileUp::origin,
            attitude.Forget<OrthogonalMap>()),
        angular_velocity_of_rigid_pile_up,
        NonRotatingPileUp::unmoving);
  }();
  if (apparent_part_rigid_motion_.empty()) {
    apparent_angular_momentum_controller_.Clear();

    for (auto& [part, actual_rigid_motion] : actual_part_rigid_motion_) {
      actual_rigid_motion =
          pile_up_actual_motion * RigidMotion<RigidPart, PileUpPrincipalAxes>(
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

  MechanicalSystem<Apparent, ApparentPileUp> apparent_system;
  for (auto const& [part, apparent_part_rigid_motion] :
       apparent_part_rigid_motion_) {
    apparent_system.AddRigidBody(
        apparent_part_rigid_motion, part->mass(), part->inertia_tensor());
  }
  auto const apparent_centre_of_mass = apparent_system.centre_of_mass();
  auto const apparent_angular_momentum = apparent_system.AngularMomentum();
  auto const apparent_inertia_tensor = apparent_system.InertiaTensor();

  // TODO(egg): Guard against the discontinuity when the order of the principal
  // axes changes.
  using SignedPileUpPrincipalAxes =
      Frame<enum class SignedPileUpPrincipalAxesTag>;
  Rotation<ApparentPileUp, SignedPileUpPrincipalAxes> const
      signed_apparent_attitude =
      apparent_inertia_tensor.Diagonalize<SignedPileUpPrincipalAxes>()
          .rotation.Inverse();
  Rotation<PileUpPrincipalAxes, NonRotatingPileUp> const actual_attitude =
      pile_up_actual_motion.orthogonal_map().AsRotation();

  Angle α = Infinity<Angle>;
  std::optional<Rotation<ApparentPileUp, PileUpPrincipalAxes>>
      apparent_attitude;

  for (auto const x : {Sign::Positive(), Sign::Negative()}) {
    for (auto const y : {Sign::Positive(), Sign::Negative()}) {
      Signature<SignedPileUpPrincipalAxes, PileUpPrincipalAxes> s(
          x, y, DeduceSignPreservingOrientation{});
      Angle const tentative_α = Abs(
          (actual_attitude * s.Forget<Rotation>() * signed_apparent_attitude)
              .RotationAngle());
      if (tentative_α < α) {
        α = tentative_α;
        apparent_attitude = s.Forget<Rotation>() * signed_apparent_attitude;
      }
    }
  }

  // The angular velocity of a rigid body with the inertia and angular momentum
  // of the apparent parts.
  auto const apparent_equivalent_angular_velocity =
      apparent_angular_momentum / apparent_inertia_tensor;

  RigidMotion<ApparentPileUp, PileUpPrincipalAxes> const
      apparent_pile_up_motion(
          RigidTransformation<ApparentPileUp, PileUpPrincipalAxes>(
              ApparentPileUp::origin,
              PileUpPrincipalAxes::origin,
              apparent_attitude->Forget<OrthogonalMap>()),
          apparent_equivalent_angular_velocity,
          ApparentPileUp::unmoving);
  RigidMotion<ApparentPileUp, NonRotatingPileUp> const
      apparent_to_pile_up_rotational_motion =
          conserve_angular_momentum
              ? pile_up_actual_motion * apparent_pile_up_motion
              : RigidMotion<ApparentPileUp, NonRotatingPileUp>::Identity();
  RigidMotion<Apparent, NonRotatingPileUp> const apparent_to_pile_up_motion =
      apparent_to_pile_up_rotational_motion *
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

  AngularVelocity<NonRotatingPileUp> const ω_actual =
      pile_up_actual_motion.Inverse().angular_velocity_of_to_frame();
      
  logger_.Append("conserveAngularMomentum", conserve_angular_momentum);
  logger_.Append("bodyFixedForces", body_fixed_forces);
  logger_.Append("inertiallyFixedForces", inertially_fixed_forces);
  logger_.Append("precalculatedTorque", precalculated_torque);
  logger_.Append("instantInitialForces", instant_initial_forces);
  logger_.Append(
      "angularVelocity",
      ω_actual,
      mathematica::ExpressIn(2 * π * Radian, quantities::si::Minute));

  std::stringstream s;
  s << "|Lap|: " << apparent_angular_momentum.Norm() << "\n"
    << "|Lac|: " << angular_momentum_.Norm() << "\n"
    << "|Lap-Lac|: "
    << (angular_momentum_ - Identity<ApparentPileUp, NonRotatingPileUp>()(
                                apparent_angular_momentum))
           .Norm()
    << "\n"
    << "|Lap|-|Lac|: "
    << angular_momentum_.Norm() - apparent_angular_momentum.Norm() << "\n"
    << u8"∡Lap, Lac: "
    << geometry::AngleBetween(angular_momentum_,
                              Identity<ApparentPileUp, NonRotatingPileUp>()(
                                  apparent_angular_momentum)) /
           quantities::si::Degree
    << u8"°\n"
    << u8"α: " << α / quantities::si::Degree << u8"°\n"
    << u8"|ωap|: "
    << apparent_equivalent_angular_velocity.Norm() /
           (2 * π * Radian / quantities::si::Minute)
    << " rpm\n"
    << u8"|ωac|: "
    << ω_actual.Norm() / (2 * π * Radian / quantities::si::Minute) << " rpm\n"
    << u8"|ωac|-|ωap|: "
    << (ω_actual.Norm() - apparent_equivalent_angular_velocity.Norm()) /
           (2 * π * Radian / quantities::si::Minute)
    << " rpm\n"
    << u8"∡ωac, ωap: "
    << geometry::AngleBetween(
           actual_attitude.Inverse()(ω_actual),
           (*apparent_attitude)(apparent_equivalent_angular_velocity)) /
           quantities::si::Degree
    << u8"°\n";
  trace = s.str();

  // This is what takes the change in moment of inertia into account.
  MakeEulerSolver(apparent_to_pile_up_rotational_motion.orthogonal_map()(
                      apparent_inertia_tensor),
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

bool PileUp::conserve_angular_momentum = true;
bool PileUp::body_fixed_forces = false;
bool PileUp::inertially_fixed_forces = false;
bool PileUp::precalculated_torque = false;
bool PileUp::instant_initial_forces = true;

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
