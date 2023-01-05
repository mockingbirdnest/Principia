#pragma once

#include "physics/ephemeris.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "absl/container/btree_set.h"
#include "absl/strings/str_cat.h"
#include "astronomy/epoch.hpp"
#include "base/jthread.hpp"
#include "base/macros.hpp"
#include "base/map_util.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "integrators/integrators.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/hermite3.hpp"
#include "physics/continuous_trajectory.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_ephemeris {

using astronomy::J2000;
using base::dynamic_cast_not_null;
using base::FindOrDie;
using base::make_not_null_unique;
using base::MakeStoppableThread;
using geometry::Barycentre;
using geometry::Displacement;
using geometry::InnerProduct;
using geometry::Position;
using geometry::R3Element;
using geometry::Sign;
using geometry::Velocity;
using integrators::EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator;
using integrators::ExplicitSecondOrderOrdinaryDifferentialEquation;
using integrators::InitialValueProblem;
using integrators::Integrator;
using integrators::methods::Fine1987RKNG34;
using numerics::Brent;
using numerics::DoublePrecision;
using numerics::Hermite3;
using quantities::Abs;
using quantities::Exponentiation;
using quantities::GravitationalParameter;
using quantities::Inverse;
using quantities::Quotient;
using quantities::Sqrt;
using quantities::Square;
using quantities::Time;
using quantities::Variation;
using quantities::si::Day;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

using namespace std::chrono_literals;

constexpr Length pre_ἐρατοσθένης_default_ephemeris_fitting_tolerance =
    1 * Milli(Metre);
constexpr Time max_time_between_checkpoints = 180 * Day;
// Below this threshold detect a collision to prevent the integrator and the
// downsampling from going postal.
constexpr double min_radius_tolerance = 0.99;

inline absl::Status CollisionDetected() {
  return absl::OutOfRangeError("Collision detected");
}

template<typename Frame>
Ephemeris<Frame>::AccuracyParameters::AccuracyParameters(
    Length const& fitting_tolerance,
    double const geopotential_tolerance)
    : fitting_tolerance_(fitting_tolerance),
      geopotential_tolerance_(geopotential_tolerance) {}

template<typename Frame>
void Ephemeris<Frame>::AccuracyParameters::WriteToMessage(
    not_null<serialization::Ephemeris::AccuracyParameters*> const message)
    const {
  fitting_tolerance_.WriteToMessage(message->mutable_fitting_tolerance());
  message->set_geopotential_tolerance(geopotential_tolerance_);
}

template<typename Frame>
typename Ephemeris<Frame>::AccuracyParameters
Ephemeris<Frame>::AccuracyParameters::ReadFromMessage(
    serialization::Ephemeris::AccuracyParameters const& message) {
  return AccuracyParameters(
      Length::ReadFromMessage(message.fitting_tolerance()),
      message.geopotential_tolerance());
}

template<typename Frame>
Ephemeris<Frame>::Ephemeris(
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies,
    std::vector<DegreesOfFreedom<Frame>> const& initial_state,
    Instant const& initial_time,
    AccuracyParameters const& accuracy_parameters,
    FixedStepParameters fixed_step_parameters)
    : accuracy_parameters_(accuracy_parameters),
      fixed_step_parameters_(std::move(fixed_step_parameters)),
      checkpointer_(
          make_not_null_unique<Checkpointer<serialization::Ephemeris>>(
              MakeCheckpointerWriter(),
              MakeCheckpointerReader())),
      reanimator_(
          [this](Instant const& desired_t_min) {
            return Reanimate(desired_t_min);
          },
          20ms) {  // 50 Hz.
  CHECK(!bodies.empty());
  CHECK_EQ(bodies.size(), initial_state.size());

  InitialValueProblem<NewtonianMotionEquation> problem;
  problem.equation = MakeMassiveBodiesNewtonianMotionEquation();

  typename NewtonianMotionEquation::SystemState& state = problem.initial_state;
  state.time = DoublePrecision<Instant>(initial_time);

  for (int i = 0; i < bodies.size(); ++i) {
    auto& body = bodies[i];
    DegreesOfFreedom<Frame> const& degrees_of_freedom = initial_state[i];

    unowned_bodies_.emplace_back(body.get());
    unowned_bodies_indices_.emplace(body.get(), i);

    auto const [it, inserted] = bodies_to_trajectories_.emplace(
        body.get(),
        std::make_unique<ContinuousTrajectory<Frame>>(
            fixed_step_parameters_.step(),
            accuracy_parameters_.fitting_tolerance_));
    CHECK(inserted);
    ContinuousTrajectory<Frame>* const trajectory = it->second.get();
    CHECK_OK(trajectory->Append(initial_time, degrees_of_freedom));

    if (body->is_oblate()) {
      geopotentials_.emplace(
          geopotentials_.cbegin(),
          dynamic_cast_not_null<OblateBody<Frame> const*>(body.get()),
          accuracy_parameters_.geopotential_tolerance_);
      // Inserting at the beginning of the vectors is O(N).
      bodies_.insert(bodies_.begin(), std::move(body));
      trajectories_.insert(trajectories_.begin(), trajectory);
      state.positions.emplace(state.positions.begin(),
                              degrees_of_freedom.position());
      state.velocities.emplace(state.velocities.begin(),
                               degrees_of_freedom.velocity());
      ++number_of_oblate_bodies_;
    } else {
      // Inserting at the end of the vectors is O(1).
      bodies_.push_back(std::move(body));
      trajectories_.push_back(trajectory);
      state.positions.emplace_back(degrees_of_freedom.position());
      state.velocities.emplace_back(degrees_of_freedom.velocity());
      ++number_of_spherical_bodies_;
    }
  }

  absl::ReaderMutexLock l(&lock_);  // For locking checks.
  instance_ = fixed_step_parameters_.integrator().NewInstance(
      problem,
      /*append_state=*/std::bind(
          &Ephemeris::AppendMassiveBodiesState, this, _1),
      fixed_step_parameters_.step());
}

template<typename Frame>
Ephemeris<Frame>::~Ephemeris() {
  reanimator_.Stop();
}

template<typename Frame>
std::vector<not_null<MassiveBody const*>> const&
Ephemeris<Frame>::bodies() const {
  return unowned_bodies_;
}

template<typename Frame>
not_null<ContinuousTrajectory<Frame> const*> Ephemeris<Frame>::trajectory(
    not_null<MassiveBody const*> body) const {
  return FindOrDie(bodies_to_trajectories_, body).get();
}

template<typename Frame>
bool Ephemeris<Frame>::empty() const {
  for (auto const& [_, trajectory] : bodies_to_trajectories_) {
    if (trajectory->empty()) {
      return true;
    }
  }
  return false;
}

template<typename Frame>
Instant Ephemeris<Frame>::t_min() const {
  absl::ReaderMutexLock l(&lock_);
  return t_min_locked();
}

template<typename Frame>
Instant Ephemeris<Frame>::t_max() const {
  absl::ReaderMutexLock l(&lock_);
  return t_max_locked();
}

template<typename Frame>
FixedStepSizeIntegrator<
    typename Ephemeris<Frame>::NewtonianMotionEquation> const&
Ephemeris<Frame>::planetary_integrator() const {
  return fixed_step_parameters_.integrator();
}

template<typename Frame>
absl::Status Ephemeris<Frame>::last_severe_integration_status() const {
  absl::ReaderMutexLock l(&lock_);
  return last_severe_integration_status_;
}

template<typename Frame>
void Ephemeris<Frame>::RequestReanimation(Instant const& desired_t_min) {
  reanimator_.Start();

  bool must_restart;
  {
    absl::MutexLock l(&lock_);

    // If the reanimator is asked to do significantly less work (as defined by
    // the time between checkpoints) than it is currently doing, interrupt it.
    // Note that this is fundamentally racy: for instance the reanimator may not
    // have picked the last input given by Put.  But it helps if the user was
    // doing a very long reanimation and wants to shorten it.
    must_restart = last_desired_t_min_.has_value() &&
                   last_desired_t_min_.value() + max_time_between_checkpoints <
                       desired_t_min;
    LOG_IF(WARNING, must_restart)
        << "Restarting reanimator because desired t_min went from "
        << last_desired_t_min_.value() << " to " << desired_t_min;
    last_desired_t_min_ = desired_t_min;
  }

  // Don't hold the lock while restarting, the reanimator needs it.
  if (must_restart) {
    reanimator_.Restart();
  }
  reanimator_.Put(desired_t_min);
}

template<typename Frame>
void Ephemeris<Frame>::WaitForReanimation(Instant const& desired_t_min) {
  auto desired_t_min_reached = [this, desired_t_min]() {
    lock_.AssertReaderHeld();
    return t_min_locked() <= desired_t_min;
  };

  absl::ReaderMutexLock l(&lock_);
  lock_.Await(absl::Condition(&desired_t_min_reached));
}

template<typename Frame>
absl::Status Ephemeris<Frame>::Prolong(Instant const& t) {
  // Short-circuit without locking.
  if (t <= t_max()) {
    return absl::OkStatus();
  }

  // Note that |t| may be before the last time that we integrated and still
  // after |t_max()|.  In this case we want to make sure that the integrator
  // makes progress.
  Instant t_final;
  Instant const instance_time = this->instance_time();
  if (t <= instance_time) {
    t_final = instance_time + fixed_step_parameters_.step();
  } else {
    t_final = t;
  }

  // Perform the integration.  Note that we may have to iterate until |t_max()|
  // actually reaches |t| because the last series may not be fully determined
  // after the first integration.
  absl::MutexLock l(&lock_);
  while (t_max_locked() < t) {
    instance_->Solve(t_final).IgnoreError();
    RETURN_IF_STOPPED;
    t_final += fixed_step_parameters_.step();
  }

  return absl::OkStatus();
}

template<typename Frame>
not_null<std::unique_ptr<typename Integrator<
    typename Ephemeris<Frame>::NewtonianMotionEquation>::Instance>>
Ephemeris<Frame>::NewInstance(
    std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
    IntrinsicAccelerations const& intrinsic_accelerations,
    FixedStepParameters const& parameters) {
  return StoppableNewInstance(trajectories, intrinsic_accelerations, parameters)
      .value();
}

template<typename Frame>
absl::StatusOr<not_null<std::unique_ptr<typename Integrator<
    typename Ephemeris<Frame>::NewtonianMotionEquation>::Instance>>>
Ephemeris<Frame>::StoppableNewInstance(
    std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
    IntrinsicAccelerations const& intrinsic_accelerations,
    FixedStepParameters const& parameters) {
  InitialValueProblem<NewtonianMotionEquation> problem;

  problem.equation.compute_acceleration =
      [this, intrinsic_accelerations](
          Instant const& t,
          std::vector<Position<Frame>> const& positions,
          std::vector<Vector<Acceleration, Frame>>& accelerations) {
    auto const error =
        ComputeGravitationalAccelerationByAllMassiveBodiesOnMasslessBodies(
            t,
            positions,
            accelerations);
    // Add the intrinsic accelerations.
    for (int i = 0; i < intrinsic_accelerations.size(); ++i) {
      auto const intrinsic_acceleration = intrinsic_accelerations[i];
      if (intrinsic_acceleration != nullptr) {
        accelerations[i] += intrinsic_acceleration(t);
      }
    }
    return error == absl::StatusCode::kOk ? absl::OkStatus() :
                    CollisionDetected();
  };

  CHECK(!trajectories.empty());
  auto const trajectory_last_time = (*trajectories.begin())->back().time;
  problem.initial_state.time = DoublePrecision<Instant>(trajectory_last_time);
  for (auto const& trajectory : trajectories) {
    auto const& [last_time, last_degrees_of_freedom] = trajectory->back();
    CHECK_EQ(last_time, trajectory_last_time);
    problem.initial_state.positions.emplace_back(
        last_degrees_of_freedom.position());
    problem.initial_state.velocities.emplace_back(
        last_degrees_of_freedom.velocity());
  }

  auto const append_state = std::bind(
      &Ephemeris::AppendMasslessBodiesStateToTrajectories, _1, trajectories);

  // The construction of the instance may evaluate the degrees of freedom of the
  // bodies.
  Prolong(trajectory_last_time + parameters.step()).IgnoreError();
  RETURN_IF_STOPPED;

  // NOTE(phl): For some reason the May 2021 version of absl wants an explicit
  // construction here.  Unsure if the bug is in absl, VS 2019, or both.
  return absl::StatusOr<not_null<std::unique_ptr<typename Integrator<
      typename Ephemeris<Frame>::NewtonianMotionEquation>::Instance>>>(
      parameters.integrator().NewInstance(
          problem, append_state, parameters.step()));
}

template<typename Frame>
absl::Status Ephemeris<Frame>::FlowWithAdaptiveStep(
    not_null<DiscreteTrajectory<Frame>*> const trajectory,
    IntrinsicAcceleration intrinsic_acceleration,
    Instant const& t,
    AdaptiveStepParameters const& parameters,
    std::int64_t const max_ephemeris_steps) {
  auto compute_acceleration = [this, &intrinsic_acceleration](
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) {
    auto const error =
        ComputeGravitationalAccelerationByAllMassiveBodiesOnMasslessBodies(
            t,
            positions,
            accelerations);
    if (intrinsic_acceleration != nullptr) {
      accelerations[0] += intrinsic_acceleration(t);
    }
    return error == absl::StatusCode::kOk ? absl::OkStatus() :
                    CollisionDetected();
  };

  return FlowODEWithAdaptiveStep<NewtonianMotionEquation>(
             std::move(compute_acceleration),
             trajectory,
             t,
             parameters,
             max_ephemeris_steps);
}

template<typename Frame>
absl::Status Ephemeris<Frame>::FlowWithAdaptiveStep(
    not_null<DiscreteTrajectory<Frame>*> trajectory,
    GeneralizedIntrinsicAcceleration intrinsic_acceleration,
    Instant const& t,
    GeneralizedAdaptiveStepParameters const& parameters,
    std::int64_t max_ephemeris_steps) {
  auto compute_acceleration =
      [this, &intrinsic_acceleration](
          Instant const& t,
          std::vector<Position<Frame>> const& positions,
          std::vector<Velocity<Frame>> const& velocities,
          std::vector<Vector<Acceleration, Frame>>& accelerations) {
        auto const error =
            ComputeGravitationalAccelerationByAllMassiveBodiesOnMasslessBodies(
                t,
                positions,
                accelerations);
        if (intrinsic_acceleration != nullptr) {
          accelerations[0] +=
              intrinsic_acceleration(t, {positions[0], velocities[0]});
        }
        return error == absl::StatusCode::kOk ? absl::OkStatus() :
                        CollisionDetected();
      };

  return FlowODEWithAdaptiveStep<GeneralizedNewtonianMotionEquation>(
             std::move(compute_acceleration),
             trajectory,
             t,
             parameters,
             max_ephemeris_steps);
}

template<typename Frame>
absl::Status Ephemeris<Frame>::FlowWithFixedStep(
    Instant const& t,
    typename Integrator<NewtonianMotionEquation>::Instance& instance) {
  if (empty() || t > t_max()) {
    Prolong(t).IgnoreError();
    RETURN_IF_STOPPED;
  }
  if (instance.time() == DoublePrecision<Instant>(t)) {
    return absl::OkStatus();
  }

  return instance.Solve(t);
}

template<typename Frame>
Vector<Acceleration, Frame>
Ephemeris<Frame>::ComputeGravitationalAccelerationOnMasslessBody(
    Position<Frame> const& position,
    Instant const& t) const {
  std::vector<Vector<Acceleration, Frame>> accelerations(1);
  ComputeGravitationalAccelerationByAllMassiveBodiesOnMasslessBodies(
      t,
      {position},
      accelerations);

  return accelerations[0];
}

template<typename Frame>
Vector<Acceleration, Frame>
Ephemeris<Frame>::ComputeGravitationalAccelerationOnMasslessBody(
    not_null<DiscreteTrajectory<Frame>*> const trajectory,
    Instant const& t) const {
  auto const it = trajectory->find(t);
  DegreesOfFreedom<Frame> const& degrees_of_freedom = it->degrees_of_freedom;
  return ComputeGravitationalAccelerationOnMasslessBody(
             degrees_of_freedom.position(), t);
}

template<typename Frame>
Vector<Acceleration, Frame>
Ephemeris<Frame>::ComputeGravitationalAccelerationOnMassiveBody(
    not_null<MassiveBody const*> const body,
    Instant const& t) const {
  bool const body_is_oblate = body->is_oblate();

  std::vector<Position<Frame>> positions;
  std::vector<Vector<Acceleration, Frame>> accelerations(bodies_.size());
  int b1 = -1;

  // Evaluate the |positions|.  Locking is necessary to be able to call the
  // "locked" method of each trajectory.
  {
    absl::ReaderMutexLock l(&lock_);
    positions.reserve(bodies_.size());
    for (int b = 0; b < bodies_.size(); ++b) {
      auto const& current_body = bodies_[b];
      auto const& current_body_trajectory = trajectories_[b];
      if (current_body.get() == body) {
        CHECK_EQ(-1, b1);
        b1 = b;
      }
      positions.push_back(current_body_trajectory->EvaluatePositionLocked(t));
    }
    CHECK_LE(0, b1);
  }

  if (body_is_oblate) {
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        /*body1_is_oblate=*/true,
        /*body2_is_oblate=*/true>(
        t,
        /*body1=*/*body, b1,
        /*bodies2=*/bodies_,
        /*b2_begin=*/0, /*b2_end=*/b1,
        positions, accelerations, geopotentials_);
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        /*body1_is_oblate=*/true,
        /*body2_is_oblate=*/true>(
        t,
        /*body1=*/*body, b1,
        /*bodies2=*/bodies_,
        /*b2_begin=*/b1 + 1, /*b2_end=*/number_of_oblate_bodies_,
        positions, accelerations, geopotentials_);
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        /*body1_is_oblate=*/true,
        /*body2_is_oblate=*/false>(
        t,
        /*body1=*/*body, b1,
        /*bodies2=*/bodies_,
        /*b2_begin=*/number_of_oblate_bodies_,
        /*b2_end=*/number_of_oblate_bodies_ + number_of_spherical_bodies_,
        positions, accelerations, geopotentials_);
  } else {
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        /*body1_is_oblate=*/false,
        /*body2_is_oblate=*/true>(
        t,
        /*body1=*/*body, b1,
        /*bodies2=*/bodies_,
        /*b2_begin=*/0, /*b2_end=*/number_of_oblate_bodies_,
        positions, accelerations, geopotentials_);
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        /*body1_is_oblate=*/false,
        /*body2_is_oblate=*/false>(
        t,
        /*body1=*/*body, b1,
        /*bodies2=*/bodies_,
        /*b2_begin=*/number_of_oblate_bodies_,
        /*b2_end=*/b1,
        positions, accelerations, geopotentials_);
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        /*body1_is_oblate=*/false,
        /*body2_is_oblate=*/false>(
        t,
        /*body1=*/*body, b1,
        /*bodies2=*/bodies_,
        /*b2_begin=*/b1 + 1,
        /*b2_end=*/number_of_oblate_bodies_ + number_of_spherical_bodies_,
        positions, accelerations, geopotentials_);
  }

  return accelerations[b1];
}

template<typename Frame>
SpecificEnergy Ephemeris<Frame>::ComputeGravitationalPotential(
    Position<Frame> const& position,
    Instant const& t) const {
  std::vector<SpecificEnergy> potentials(1);
  ComputeGravitationalPotentialsOfAllMassiveBodies(t, {position}, potentials);

  return potentials[0];
}

template<typename Frame>
void Ephemeris<Frame>::ComputeApsides(not_null<MassiveBody const*> const body1,
                                      not_null<MassiveBody const*> const body2,
                                      DiscreteTrajectory<Frame>& apoapsides1,
                                      DiscreteTrajectory<Frame>& periapsides1,
                                      DiscreteTrajectory<Frame>& apoapsides2,
                                      DiscreteTrajectory<Frame>& periapsides2) {
  absl::ReaderMutexLock l(&lock_);

  not_null<ContinuousTrajectory<Frame> const*> const body1_trajectory =
      trajectory(body1);
  not_null<ContinuousTrajectory<Frame> const*> const body2_trajectory =
      trajectory(body2);

  // Computes the derivative of the squared distance between |body1| and |body2|
  // at time |t|.
  auto const evaluate_square_distance_derivative =
      [body1_trajectory, body2_trajectory](
          Instant const& t) -> Variation<Square<Length>> {
    DegreesOfFreedom<Frame> const body1_degrees_of_freedom =
        body1_trajectory->EvaluateDegreesOfFreedomLocked(t);
    DegreesOfFreedom<Frame> const body2_degrees_of_freedom =
        body2_trajectory->EvaluateDegreesOfFreedomLocked(t);
    RelativeDegreesOfFreedom<Frame> const relative =
        body1_degrees_of_freedom - body2_degrees_of_freedom;
    return 2.0 * InnerProduct(relative.displacement(), relative.velocity());
  };

  std::optional<Instant> previous_time;
  std::optional<Variation<Square<Length>>> previous_squared_distance_derivative;

  for (Instant time = t_min_locked();
       time <= t_max_locked();
       time += fixed_step_parameters_.step()) {
    Variation<Square<Length>> const squared_distance_derivative =
        evaluate_square_distance_derivative(time);
    if (previous_squared_distance_derivative &&
        Sign(squared_distance_derivative) !=
            Sign(*previous_squared_distance_derivative)) {
      CHECK(previous_time);

      // The derivative of |squared_distance| changed sign.  Find its zero by
      // Brent's method, this is the time of the apsis.  Then compute the apsis
      // and append it to one of the output trajectories.
      Instant const apsis_time = Brent(evaluate_square_distance_derivative,
                                       *previous_time,
                                       time);
      DegreesOfFreedom<Frame> const apsis1_degrees_of_freedom =
          body1_trajectory->EvaluateDegreesOfFreedomLocked(apsis_time);
      DegreesOfFreedom<Frame> const apsis2_degrees_of_freedom =
          body2_trajectory->EvaluateDegreesOfFreedomLocked(apsis_time);
      if (Sign(squared_distance_derivative).is_negative()) {
        apoapsides1.Append(apsis_time, apsis1_degrees_of_freedom).IgnoreError();
        apoapsides2.Append(apsis_time, apsis2_degrees_of_freedom).IgnoreError();
      } else {
        periapsides1.Append(apsis_time,
                            apsis1_degrees_of_freedom).IgnoreError();
        periapsides2.Append(apsis_time,
                            apsis2_degrees_of_freedom).IgnoreError();
      }
    }

    previous_time = time;
    previous_squared_distance_derivative = squared_distance_derivative;
  }
}

template<typename Frame>
int Ephemeris<Frame>::serialization_index_for_body(
    not_null<MassiveBody const*> const body) const {
  return FindOrDie(unowned_bodies_indices_, body);
}

template<typename Frame>
not_null<MassiveBody const*> Ephemeris<Frame>::body_for_serialization_index(
    int const serialization_index) const {
  return unowned_bodies_[serialization_index];
}

template<typename Frame>
void Ephemeris<Frame>::WriteToMessage(
    not_null<serialization::Ephemeris*> const message) const {
  LOG(INFO) << __FUNCTION__;
  absl::ReaderMutexLock l(&lock_);

  // Make sure that a checkpoint exists, otherwise we would not serialize some
  // parts of the state.
  WriteToCheckpointIfNeeded(instance_->time().value);
  checkpointer_->WriteToMessage(message->mutable_checkpoint());

  // The bodies are serialized in the order in which they were given at
  // construction.
  for (auto const& unowned_body : unowned_bodies_) {
    unowned_body->WriteToMessage(message->add_body());
  }
  // The trajectories are serialized in the order resulting from the separation
  // between oblate and spherical bodies.
  for (auto const& trajectory : trajectories_) {
    trajectory->WriteToMessage(message->add_trajectory());
  }
  fixed_step_parameters_.WriteToMessage(
      message->mutable_fixed_step_parameters());
  accuracy_parameters_.WriteToMessage(
      message->mutable_accuracy_parameters());
  LOG(INFO) << NAMED(message->SpaceUsed());
  LOG(INFO) << NAMED(message->ByteSize());
}

template<typename Frame>
template<typename, typename>
not_null<std::unique_ptr<Ephemeris<Frame>>> Ephemeris<Frame>::ReadFromMessage(
    Instant const& desired_t_min,
    serialization::Ephemeris const& message) {
  bool const is_pre_ἐρατοσθένης = !message.has_accuracy_parameters();
  bool const is_pre_fatou = !message.has_checkpoint_time();
  bool const is_pre_grassmann = message.checkpoint_size() == 0;
  LOG_IF(WARNING, is_pre_grassmann)
      << "Reading pre-"
      << (is_pre_ἐρατοσθένης ? "Ἐρατοσθένης"
          : is_pre_fatou     ? "Fatou"
                             : "Grassmann") << " Ephemeris";

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  for (auto const& body : message.body()) {
    bodies.push_back(MassiveBody::ReadFromMessage(body));
  }

  AccuracyParameters accuracy_parameters(
      pre_ἐρατοσθένης_default_ephemeris_fitting_tolerance,
      /*geopotential_tolerance=*/0);
  if (!is_pre_ἐρατοσθένης) {
    accuracy_parameters =
        AccuracyParameters::ReadFromMessage(message.accuracy_parameters());
  }
  FixedStepParameters const fixed_step_parameters =
      FixedStepParameters::ReadFromMessage(message.fixed_step_parameters());

  // Dummy initial state and time.  We'll overwrite them later.
  std::vector<DegreesOfFreedom<Frame>> const initial_state(
      bodies.size(),
      DegreesOfFreedom<Frame>(Frame::origin, Frame::unmoving));
  Instant const initial_time;
  auto ephemeris = make_not_null_unique<Ephemeris<Frame>>(
                       std::move(bodies),
                       initial_state,
                       initial_time,
                       accuracy_parameters,
                       fixed_step_parameters);

  int index = 0;
  ephemeris->bodies_to_trajectories_.clear();
  ephemeris->trajectories_.clear();
  for (auto const& trajectory : message.trajectory()) {
    not_null<MassiveBody const*> const body = ephemeris->bodies_[index].get();
    not_null<std::unique_ptr<ContinuousTrajectory<Frame>>>
        deserialized_trajectory = ContinuousTrajectory<Frame>::ReadFromMessage(
            desired_t_min, trajectory);
    ephemeris->trajectories_.push_back(deserialized_trajectory.get());
    ephemeris->bodies_to_trajectories_.emplace(
        body, std::move(deserialized_trajectory));
    ++index;
  }
  CHECK_LT(0, index) << "Empty ephemeris";

  if (is_pre_grassmann) {
    serialization::Ephemeris serialized_ephemeris;
    auto* const checkpoint = serialized_ephemeris.add_checkpoint();
    if (is_pre_fatou) {
      *checkpoint->mutable_time() =
          message.instance().current_state().time().value().point();
    } else {
      *checkpoint->mutable_time() = message.checkpoint_time();
    }
    *checkpoint->mutable_instance() = message.instance();
    ephemeris->checkpointer_ =
        Checkpointer<serialization::Ephemeris>::ReadFromMessage(
            ephemeris->MakeCheckpointerWriter(),
            ephemeris->MakeCheckpointerReader(),
            serialized_ephemeris.checkpoint());
  } else {
    ephemeris->checkpointer_ =
        Checkpointer<serialization::Ephemeris>::ReadFromMessage(
            ephemeris->MakeCheckpointerWriter(),
            ephemeris->MakeCheckpointerReader(),
            message.checkpoint());
  }

  // The checkpoint at or before |desired_t_min| will result in a |t_min()|
  // which is at |desired_t_min| (if the checkpoint was taken with
  // |last_points_.size() == 1|) or before (if the checkpoint was taken with
  // |last_points_.size() > 1|).
  ephemeris->oldest_reanimated_checkpoint_ =
      ephemeris->checkpointer_->checkpoint_at_or_before(desired_t_min);
  if (ephemeris->oldest_reanimated_checkpoint_ == InfinitePast) {
    // In the pre-Grassmann compatibility case the (only) checkpoint may be
    // after |desired_t_min|.  This also happens with old saves that are
    // rewritten post-Grassmann.
    CHECK_LE(ephemeris->t_min(), desired_t_min);
  } else {
    LOG(INFO) << "Restoring to checkpoint at "
              << ephemeris->oldest_reanimated_checkpoint_;
    CHECK_OK(ephemeris->checkpointer_->ReadFromCheckpointAt(
        ephemeris->oldest_reanimated_checkpoint_));
  }

  // The ephemeris will need to be prolonged and reanimated as needed when
  // deserializing the plugin.
  return ephemeris;
}

template<typename Frame>
Ephemeris<Frame>::Ephemeris(
    FixedStepSizeIntegrator<
        typename Ephemeris<Frame>::NewtonianMotionEquation> const& integrator)
    : accuracy_parameters_(pre_ἐρατοσθένης_default_ephemeris_fitting_tolerance,
                           /*geopotential_tolerance=*/0),
      fixed_step_parameters_(integrator, 1 * Second),
      checkpointer_(
          make_not_null_unique<Checkpointer<serialization::Ephemeris>>(
              /*reader=*/nullptr, /*writer=*/nullptr)),
      reanimator_(/*action=*/nullptr, 0ms) {}

template<typename Frame>
void Ephemeris<Frame>::WriteToCheckpointIfNeeded(Instant const& time) const {
  if constexpr (base::is_serializable_v<Frame>) {
    lock_.AssertReaderHeld();
    if (checkpointer_->WriteToCheckpointIfNeeded(
            time, max_time_between_checkpoints)) {
      for (auto const& trajectory : trajectories_) {
        trajectory->WriteToCheckpoint(time);
      }
    }
  }
}

template<typename Frame>
Checkpointer<serialization::Ephemeris>::Writer
Ephemeris<Frame>::MakeCheckpointerWriter() {
  if constexpr (base::is_serializable_v<Frame>) {
    return [this](
               not_null<serialization::Ephemeris::Checkpoint*> const message) {
      lock_.AssertReaderHeld();
      instance_->WriteToMessage(message->mutable_instance());
    };
  } else {
    return nullptr;
  }
}

template<typename Frame>
Checkpointer<serialization::Ephemeris>::Reader
Ephemeris<Frame>::MakeCheckpointerReader() {
  if constexpr (base::is_serializable_v<Frame>) {
    return [this](serialization::Ephemeris::Checkpoint const& message) {
      absl::MutexLock l(&lock_);
      instance_ = FixedStepSizeIntegrator<NewtonianMotionEquation>::Instance::
          ReadFromMessage(
              message.instance(),
              MakeMassiveBodiesNewtonianMotionEquation(),
              /*append_state=*/
              std::bind(&Ephemeris::AppendMassiveBodiesState, this, _1));
      return absl::OkStatus();
    };
  } else {
    return nullptr;
  }
}

template<typename Frame>
absl::Status Ephemeris<Frame>::Reanimate(Instant const desired_t_min) {
  absl::btree_set<Instant> checkpoints;
  {
    absl::ReaderMutexLock l(&lock_);

    // It is very important that |oldest_reanimated_checkpoint_| be only read by
    // the |reanimator_| thread.  If the caller was trying to determine the set
    // of checkpoints to reanimate it might race with a reanimation already in
    // flight and result in the same checkpoint reanimated multiple times, which
    // is a no-no.
    Instant const oldest_checkpoint_to_reanimate =
        checkpointer_->checkpoint_at_or_before(desired_t_min);
    checkpoints = checkpointer_->all_checkpoints_between(
        oldest_checkpoint_to_reanimate, oldest_reanimated_checkpoint_);
  }

  // This loop integrates all the segments defined by the checkpoints, going
  // backwards in time.  The last checkpoint is not restored, it just serves as
  // a limit.
  std::optional<Instant> following_checkpoint;
  for (auto it = checkpoints.crbegin(); it != checkpoints.crend(); ++it) {
    Instant const& checkpoint = *it;
    if (following_checkpoint.has_value()) {
      RETURN_IF_ERROR(checkpointer_->ReadFromCheckpointAt(
          checkpoint,
          [this,
           t_final = following_checkpoint.value(),
           t_initial = checkpoint](
              serialization::Ephemeris::Checkpoint const& message) {
            if constexpr (base::is_serializable_v<Frame>) {
              return ReanimateOneCheckpoint(message, t_initial, t_final);
            } else {
              return absl::UnknownError(
                  "No reanimation for non-serializable frames");
            }
          }));
    }
    following_checkpoint = checkpoint;
  }
  return absl::OkStatus();
}

template<typename Frame>
absl::Status Ephemeris<Frame>::ReanimateOneCheckpoint(
    serialization::Ephemeris::Checkpoint const& message,
    Instant const& t_initial,
    Instant const& t_final) {
  LOG(INFO) << "Reanimating segment from " << t_initial << " to " << t_final;

  // Create new trajectories and initialize them from the checkpoint at
  // t_initial.
  std::vector<not_null<std::unique_ptr<ContinuousTrajectory<Frame>>>>
      trajectories;
  for (int i = 0; i < trajectories_.size(); ++i) {
    trajectories.emplace_back(std::make_unique<ContinuousTrajectory<Frame>>(
        fixed_step_parameters_.step(),
        accuracy_parameters_.fitting_tolerance_));

    // This statement is subtle: it restores the checkpoints of the trajectories
    // of this ephemeris, but thanks to the newly-created reader, it restores
    // them into the local trajectories.
    CHECK_OK(trajectories_[i]->ReadFromCheckpointAt(
        t_initial, trajectories[i]->MakeCheckpointerReader()));
  }

  // Reconstruct the integrator instance from the current checkpoint.
  auto append_massive_bodies_state =
      [&trajectories](
          typename NewtonianMotionEquation::SystemState const& state) {
        AppendMassiveBodiesStateToTrajectories(state, trajectories);
      };
  auto const instance = FixedStepSizeIntegrator<NewtonianMotionEquation>::
      Instance::ReadFromMessage(message.instance(),
                                MakeMassiveBodiesNewtonianMotionEquation(),
                                append_massive_bodies_state);

  // Do the integration.  After this step the t_max() of the trajectories may
  // be before t_final because there may be last_points_ that haven't been put
  // in a series.  Don't proceed in case of error, we would run into a gap when
  // trying to stitch the trajectories.
  RETURN_IF_ERROR(instance->Solve(t_final));

  // Stitch the local trajectories to the ones in this ephemeris and record that
  // we will not reanimate this checkpoint again.
  {
    absl::MutexLock l(&lock_);
    for (int i = 0; i < trajectories_.size(); ++i) {
      trajectories_[i]->Prepend(std::move(*trajectories[i]));
    }
    oldest_reanimated_checkpoint_ = t_initial;
  }

  return absl::OkStatus();
}

template<typename Frame>
void Ephemeris<Frame>::AppendMassiveBodiesState(
    typename NewtonianMotionEquation::SystemState const& state) {
  lock_.AssertHeld();

  // Extend the trajectories.
  auto const statuses = AppendMassiveBodiesStateToTrajectories(state,
                                                               trajectories_);

  // Handle the apocalypse.
  for (int i = 0; i < statuses.size(); ++i) {
    auto const& status = statuses[i];
    if (!status.ok()) {
      last_severe_integration_status_ =
          absl::Status(status.code(),
                       absl::StrCat("Error extending trajectory for ",
                                    bodies_[i]->name(), ". ",
                                    status.message()));
      LOG(ERROR) << "New Apocalypse: " << last_severe_integration_status_;
    }
  }

  // Note that the checkpoint is written systematically after inserting the
  // first point of the trajectories.
  WriteToCheckpointIfNeeded(state.time.value);
}

template<typename Frame>
template<typename ContinuousTrajectoryPtr>
std::vector<absl::Status>
Ephemeris<Frame>::AppendMassiveBodiesStateToTrajectories(
    typename NewtonianMotionEquation::SystemState const& state,
    std::vector<not_null<ContinuousTrajectoryPtr>> const& trajectories) {
  std::vector<absl::Status> statuses;
  statuses.reserve(trajectories.size());
  Instant const time = state.time.value;
  int index = 0;
  for (auto& trajectory : trajectories) {
    statuses.push_back(trajectory->Append(
        time,
        DegreesOfFreedom<Frame>(state.positions[index].value,
                                state.velocities[index].value)));
    ++index;
  }
  return statuses;
}

template<typename Frame>
void Ephemeris<Frame>::AppendMasslessBodiesStateToTrajectories(
    typename NewtonianMotionEquation::SystemState const& state,
    std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories) {
  Instant const time = state.time.value;
  int index = 0;
  for (auto& trajectory : trajectories) {
    trajectory->Append(
        time,
        DegreesOfFreedom<Frame>(state.positions[index].value,
                                state.velocities[index].value)).IgnoreError();
    ++index;
  }
}

template<typename Frame>
typename Ephemeris<Frame>::NewtonianMotionEquation
Ephemeris<Frame>::MakeMassiveBodiesNewtonianMotionEquation() {
  NewtonianMotionEquation equation;
  equation.compute_acceleration =
      [this](Instant const& t,
             std::vector<Position<Frame>> const& positions,
             std::vector<Vector<Acceleration, Frame>>& accelerations) {
        return ComputeGravitationalAccelerationBetweenAllMassiveBodies(
                   t,
                   positions,
                   accelerations);
      };
  return equation;
}

template<typename Frame>
Instant Ephemeris<Frame>::instance_time() const {
  absl::ReaderMutexLock l(&lock_);
  return instance_->time().value;
}

template<typename Frame>
Instant Ephemeris<Frame>::t_min_locked() const {
  lock_.AssertReaderHeld();
  Instant t_min = bodies_to_trajectories_.begin()->second->t_min();
  for (auto const& [_, trajectory] : bodies_to_trajectories_) {
    t_min = std::max(t_min, trajectory->t_min_locked());
  }
  return t_min;
}

template<typename Frame>
Instant Ephemeris<Frame>::t_max_locked() const {
  lock_.AssertReaderHeld();
  Instant t_max = bodies_to_trajectories_.begin()->second->t_max();
  for (auto const& [_, trajectory] : bodies_to_trajectories_) {
    t_max = std::min(t_max, trajectory->t_max_locked());
  }
  return t_max;
}

template<typename Frame>
template<bool body1_is_oblate,
         bool body2_is_oblate,
         typename MassiveBodyConstPtr>
void Ephemeris<Frame>::
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies(
        Instant const& t,
        MassiveBody const& body1,
        std::size_t const b1,
        std::vector<not_null<MassiveBodyConstPtr>> const& bodies2,
        std::size_t const b2_begin,
        std::size_t const b2_end,
        std::vector<Position<Frame>> const& positions,
        std::vector<Vector<Acceleration, Frame>>& accelerations,
        std::vector<Geopotential<Frame>> const& geopotentials) {
  Position<Frame> const& position_of_b1 = positions[b1];
  Vector<Acceleration, Frame>& acceleration_on_b1 = accelerations[b1];
  GravitationalParameter const& μ1 = body1.gravitational_parameter();
  for (std::size_t b2 = b2_begin; b2 < b2_end; ++b2) {
    Vector<Acceleration, Frame>& acceleration_on_b2 = accelerations[b2];
    MassiveBody const& body2 = *bodies2[b2];
    GravitationalParameter const& μ2 = body2.gravitational_parameter();

    // A vector from the center of |b2| to the center of |b1|.
    Displacement<Frame> const Δq = position_of_b1 - positions[b2];

    Square<Length> const Δq² = Δq.Norm²();
    Length const Δq_norm = Sqrt(Δq²);
    Exponentiation<Length, -3> const one_over_Δq³ = Δq_norm / (Δq² * Δq²);

    auto const μ1_over_Δq³ = μ1 * one_over_Δq³;
    acceleration_on_b2 += Δq * μ1_over_Δq³;

    // [New87], Lex. III. Actioni contrariam semper & æqualem esse reactionem:
    // sive corporum duorum actiones in se mutuo semper esse æquales &
    // in partes contrarias dirigi.
    auto const μ2_over_Δq³ = μ2 * one_over_Δq³;
    acceleration_on_b1 -= Δq * μ2_over_Δq³;

    if (body1_is_oblate || body2_is_oblate) {
      if (body1_is_oblate) {
        Vector<Quotient<Acceleration,
                        GravitationalParameter>, Frame> const
            spherical_harmonics_effect =
                geopotentials[b1].GeneralSphericalHarmonicsAcceleration(
                    t,
                    -Δq,
                    Δq_norm,
                    Δq²,
                    one_over_Δq³);
        acceleration_on_b1 -= μ2 * spherical_harmonics_effect;
        acceleration_on_b2 += μ1 * spherical_harmonics_effect;
      }
      if (body2_is_oblate) {
        Vector<Quotient<Acceleration,
                        GravitationalParameter>, Frame> const
            degree_2_zonal_effect2 =
                geopotentials[b2].GeneralSphericalHarmonicsAcceleration(
                    t,
                    Δq,
                    Δq_norm,
                    Δq²,
                    one_over_Δq³);
        acceleration_on_b1 += μ2 * degree_2_zonal_effect2;
        acceleration_on_b2 -= μ1 * degree_2_zonal_effect2;
      }
    }
  }
}

template<typename Frame>
template<bool body1_is_oblate>
std::underlying_type_t<absl::StatusCode>
Ephemeris<Frame>::ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies(
    Instant const& t,
    MassiveBody const& body1,
    std::size_t const b1,
    std::vector<Position<Frame>> const& positions,
    std::vector<Vector<Acceleration, Frame>>& accelerations) const {
  lock_.AssertReaderHeld();
  GravitationalParameter const& μ1 = body1.gravitational_parameter();
  auto const& trajectory1 = *trajectories_[b1];
  Position<Frame> const position1 = trajectory1.EvaluatePositionLocked(t);
  Length const body1_collision_radius =
      min_radius_tolerance * body1.min_radius();
  // TODO(phl): Use std::to_underlying when we have C++23.
  auto error = static_cast<std::underlying_type_t<absl::StatusCode>>(
      absl::StatusCode::kOk);

  for (std::size_t b2 = 0; b2 < positions.size(); ++b2) {
    // A vector from the center of |b2| to the center of |b1|.
    Displacement<Frame> const Δq = position1 - positions[b2];

    Square<Length> const Δq² = Δq.Norm²();
    Length const Δq_norm = Sqrt(Δq²);
    error |= Δq_norm > body1_collision_radius
                 ? static_cast<std::underlying_type_t<absl::StatusCode>>(
                       absl::StatusCode::kOk)
                 : static_cast<std::underlying_type_t<absl::StatusCode>>(
                       absl::StatusCode::kOutOfRange);

    Exponentiation<Length, -3> const one_over_Δq³ = Δq_norm / (Δq² * Δq²);

    auto const μ1_over_Δq³ = μ1 * one_over_Δq³;
    accelerations[b2] += Δq * μ1_over_Δq³;

    if (body1_is_oblate) {
      Vector<Quotient<Acceleration,
                      GravitationalParameter>, Frame> const
          spherical_harmonics_effect =
              geopotentials_[b1].GeneralSphericalHarmonicsAcceleration(
                  t,
                  -Δq,
                  Δq_norm,
                  Δq²,
                  one_over_Δq³);
      accelerations[b2] += μ1 * spherical_harmonics_effect;
    }
  }
  return error;
}

template<typename Frame>
template<bool body1_is_oblate>
void Ephemeris<Frame>::ComputeGravitationalPotentialsOfMassiveBody(
    Instant const& t,
    MassiveBody const& body1,
    std::size_t b1,
    std::vector<Position<Frame>> const& positions,
    std::vector<SpecificEnergy>& potentials) const {
  lock_.AssertReaderHeld();
  GravitationalParameter const& μ1 = body1.gravitational_parameter();
  auto const& trajectory1 = *trajectories_[b1];
  Position<Frame> const position1 = trajectory1.EvaluatePositionLocked(t);

  for (std::size_t b2 = 0; b2 < positions.size(); ++b2) {
    // A vector from the center of |b2| to the center of |b1|.
    Displacement<Frame> const Δq = position1 - positions[b2];

    Square<Length> const Δq² = Δq.Norm²();
    Length const Δq_norm = Sqrt(Δq²);
    Inverse<Length> const one_over_Δq_norm = 1 / Δq_norm;
    potentials[b2] -= μ1 * one_over_Δq_norm;

    if (body1_is_oblate) {
      Exponentiation<Length, -3> const one_over_Δq³ =
          one_over_Δq_norm * one_over_Δq_norm * one_over_Δq_norm;

      auto const μ1_over_Δq³ = μ1 * one_over_Δq³;
      Quotient<SpecificEnergy, GravitationalParameter> const
          spherical_harmonics_effect =
              geopotentials_[b1].GeneralSphericalHarmonicsPotential(
                  t,
                  -Δq,
                  Δq_norm,
                  Δq²,
                  one_over_Δq³);
      potentials[b2] += μ1 * spherical_harmonics_effect;
    }
  }
}

template<typename Frame>
absl::Status
Ephemeris<Frame>::ComputeGravitationalAccelerationBetweenAllMassiveBodies(
    Instant const& t,
    std::vector<Position<Frame>> const& positions,
    std::vector<Vector<Acceleration, Frame>>& accelerations) const {
  RETURN_IF_STOPPED;

  accelerations.assign(accelerations.size(), Vector<Acceleration, Frame>());

  for (std::size_t b1 = 0; b1 < number_of_oblate_bodies_; ++b1) {
    MassiveBody const& body1 = *bodies_[b1];
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        /*body1_is_oblate=*/true,
        /*body2_is_oblate=*/true>(
        t,
        body1, b1,
        /*bodies2=*/bodies_,
        /*b2_begin=*/b1 + 1,
        /*b2_end=*/number_of_oblate_bodies_,
        positions, accelerations, geopotentials_);
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        /*body1_is_oblate=*/true,
        /*body2_is_oblate=*/false>(
        t,
        body1, b1,
        /*bodies2=*/bodies_,
        /*b2_begin=*/number_of_oblate_bodies_,
        /*b2_end=*/number_of_oblate_bodies_ + number_of_spherical_bodies_,
        positions, accelerations, geopotentials_);
  }
  for (std::size_t b1 = number_of_oblate_bodies_;
       b1 < number_of_oblate_bodies_ +
            number_of_spherical_bodies_;
       ++b1) {
    MassiveBody const& body1 = *bodies_[b1];
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        /*body1_is_oblate=*/false,
        /*body2_is_oblate=*/false>(
        t,
        body1, b1,
        /*bodies2=*/bodies_,
        /*b2_begin=*/b1 + 1,
        /*b2_end=*/number_of_oblate_bodies_ + number_of_spherical_bodies_,
        positions, accelerations, geopotentials_);
  }

  return absl::OkStatus();
}

template<typename Frame>
absl::StatusCode
Ephemeris<Frame>::
ComputeGravitationalAccelerationByAllMassiveBodiesOnMasslessBodies(
    Instant const& t,
    std::vector<Position<Frame>> const& positions,
    std::vector<Vector<Acceleration, Frame>>& accelerations) const {
  CHECK_EQ(positions.size(), accelerations.size());
  accelerations.assign(accelerations.size(), Vector<Acceleration, Frame>());
  // TODO(phl): Use std::to_underlying when we have C++23.
  auto error = static_cast<std::underlying_type_t<absl::StatusCode>>(
      absl::StatusCode::kOk);

  // Locking ensures that we see a consistent state of all the trajectories.
  absl::ReaderMutexLock l(&lock_);
  for (std::size_t b1 = 0; b1 < number_of_oblate_bodies_; ++b1) {
    MassiveBody const& body1 = *bodies_[b1];
    error |= ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies<
                 /*body1_is_oblate=*/true>(
                 t,
                 body1, b1,
                 positions,
                 accelerations);
  }
  for (std::size_t b1 = number_of_oblate_bodies_;
       b1 < number_of_oblate_bodies_ +
            number_of_spherical_bodies_;
       ++b1) {
    MassiveBody const& body1 = *bodies_[b1];
    error |= ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies<
                 /*body1_is_oblate=*/false>(
                 t,
                 body1, b1,
                 positions,
                 accelerations);
  }
  return static_cast<absl::StatusCode>(error);
}

template<typename Frame>
void Ephemeris<Frame>::ComputeGravitationalPotentialsOfAllMassiveBodies(
    Instant const& t,
    std::vector<Position<Frame>> const& positions,
    std::vector<SpecificEnergy>& potentials) const {
  CHECK_EQ(positions.size(), potentials.size());
  potentials.assign(potentials.size(), SpecificEnergy());

  // Locking ensures that we see a consistent state of all the trajectories.
  absl::ReaderMutexLock l(&lock_);
  for (std::size_t b1 = 0; b1 < number_of_oblate_bodies_; ++b1) {
    MassiveBody const& body1 = *bodies_[b1];
    ComputeGravitationalPotentialsOfMassiveBody</*body1_is_oblate=*/true>(
        t,
        body1, b1,
        positions,
        potentials);
  }
  for (std::size_t b1 = number_of_oblate_bodies_;
       b1 < number_of_oblate_bodies_ +
            number_of_spherical_bodies_;
       ++b1) {
    MassiveBody const& body1 = *bodies_[b1];
    ComputeGravitationalPotentialsOfMassiveBody</*body1_is_oblate=*/false>(
        t,
        body1, b1,
        positions,
        potentials);
  }
}


template<typename Frame>
template<typename ODE>
absl::Status Ephemeris<Frame>::FlowODEWithAdaptiveStep(
    typename ODE::RightHandSideComputation compute_acceleration,
    not_null<DiscreteTrajectory<Frame>*> trajectory,
    Instant const& t,
    physics::AdaptiveStepParameters<ODE> const& parameters,
    std::int64_t max_ephemeris_steps) {
  auto const& [trajectory_last_time,
               trajectory_last_degrees_of_freedom] = trajectory->back();
  if (trajectory_last_time == t) {
    return absl::OkStatus();
  }

  std::vector<not_null<DiscreteTrajectory<Frame>*>> const trajectories =
      {trajectory};
  // The |min| is here to prevent us from spending too much time computing the
  // ephemeris.  The |max| is here to ensure that we always try to integrate
  // forward.  We use |last_state_.time.value| because this is always finite,
  // contrary to |t_max()|, which is -∞ when |empty()|.
  Instant const t_final =
      std::min(std::max(instance_time() +
                            max_ephemeris_steps * fixed_step_parameters_.step(),
                        trajectory_last_time + fixed_step_parameters_.step()),
               t);
  Prolong(t_final).IgnoreError();
  RETURN_IF_STOPPED;

  InitialValueProblem<ODE> problem;
  problem.equation.compute_acceleration = std::move(compute_acceleration);

  problem.initial_state = {trajectory_last_time,
                           {trajectory_last_degrees_of_freedom.position()},
                           {trajectory_last_degrees_of_freedom.velocity()}};

  typename AdaptiveStepSizeIntegrator<ODE>::Parameters const
      integrator_parameters(
          /*first_time_step=*/t_final - problem.initial_state.time.value,
          /*safety_factor=*/0.9,
          parameters.max_steps(),
          /*last_step_is_exact=*/true);
  CHECK_GT(integrator_parameters.first_step, 0 * Second)
      << "Flow back to the future: " << t_final
      << " <= " << problem.initial_state.time.value;
  auto const tolerance_to_error_ratio =
      std::bind(&Ephemeris<Frame>::ToleranceToErrorRatio,
                parameters.length_integration_tolerance(),
                parameters.speed_integration_tolerance(),
                _1, _2, _3);

  typename AdaptiveStepSizeIntegrator<ODE>::AppendState append_state =
      std::bind(&Ephemeris::AppendMasslessBodiesStateToTrajectories,
                _1,
                std::cref(trajectories));
  auto const instance =
      parameters.integrator().NewInstance(problem,
                                          append_state,
                                          tolerance_to_error_ratio,
                                          integrator_parameters);
  auto status = instance->Solve(t_final);

  // We probably don't care if the vessel gets too close to the singularity, as
  // we only use this integrator for the future.  So we swallow the error.  Note
  // that a collision in the prediction or the flight plan (for which this path
  // is used) should not cause the vessel to be deleted.
  if (absl::IsOutOfRange(status)) {
    status = absl::OkStatus();
  }

  // TODO(egg): when we have events in trajectories, we should add a singularity
  // event at the end if the outcome indicates a singularity
  // (|VanishingStepSize|).  We should not have an event on the trajectory if
  // |ReachedMaximalStepCount|, since that is not a physical property, but
  // rather a self-imposed constraint.
  if (!status.ok() || t_final == t) {
    return status;
  } else {
    return absl::DeadlineExceededError("Couldn't reach " + DebugString(t) +
                                       ", stopping at " + DebugString(t_final));
  }
}

template<typename Frame>
double Ephemeris<Frame>::ToleranceToErrorRatio(
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance,
    Time const& current_step_size,
    typename NewtonianMotionEquation::SystemState const& /*state*/,
    typename NewtonianMotionEquation::SystemStateError const& error) {
  Length max_length_error;
  Speed max_speed_error;
  for (auto const& position_error : error.position_error) {
    max_length_error = std::max(max_length_error,
                                position_error.Norm());
  }
  for (auto const& velocity_error : error.velocity_error) {
    max_speed_error = std::max(max_speed_error,
                               velocity_error.Norm());
  }
  return std::min(length_integration_tolerance / max_length_error,
                  speed_integration_tolerance / max_speed_error);
}

template<typename Frame>
typename Ephemeris<Frame>::IntrinsicAccelerations const
    Ephemeris<Frame>::NoIntrinsicAccelerations;

}  // namespace internal_ephemeris
}  // namespace physics
}  // namespace principia
