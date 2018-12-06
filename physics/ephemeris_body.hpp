
#pragma once

#include "physics/ephemeris.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <optional>
#include <set>
#include <vector>

#include "astronomy/epoch.hpp"
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
using base::Error;
using base::FindOrDie;
using base::make_not_null_unique;
using geometry::Barycentre;
using geometry::Displacement;
using geometry::InnerProduct;
using geometry::Position;
using geometry::R3Element;
using geometry::Sign;
using geometry::Velocity;
using integrators::EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator;
using integrators::ExplicitSecondOrderOrdinaryDifferentialEquation;
using integrators::Integrator;
using integrators::IntegrationProblem;
using integrators::methods::Fine1987RKNG34;
using numerics::Bisect;
using numerics::DoublePrecision;
using numerics::Hermite3;
using quantities::Abs;
using quantities::Exponentiation;
using quantities::GravitationalParameter;
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

constexpr Length pre_ἐρατοσθένης_default_ephemeris_fitting_tolerance =
    1 * Milli(Metre);
constexpr Time max_time_between_checkpoints = 180 * Day;

template<typename Frame>
Ephemeris<Frame>::AccuracyParameters::AccuracyParameters(
    Length const& fitting_tolerance)
    : fitting_tolerance_(fitting_tolerance) {}

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
Ephemeris<Frame>::AdaptiveStepParameters::AdaptiveStepParameters(
    AdaptiveStepSizeIntegrator<NewtonianMotionEquation> const& integrator,
    std::int64_t const max_steps,
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance)
    : integrator_(&integrator),
      max_steps_(max_steps),
      length_integration_tolerance_(length_integration_tolerance),
      speed_integration_tolerance_(speed_integration_tolerance) {
  CHECK_LT(0, max_steps_);
  CHECK_LT(Length(), length_integration_tolerance_);
  CHECK_LT(Speed(), speed_integration_tolerance_);
}

template<typename Frame>
AdaptiveStepSizeIntegrator<
    typename Ephemeris<Frame>::NewtonianMotionEquation> const&
Ephemeris<Frame>::AdaptiveStepParameters::integrator() const {
  return *integrator_;
}

template<typename Frame>
std::int64_t Ephemeris<Frame>::AdaptiveStepParameters::max_steps() const {
  return max_steps_;
}

template<typename Frame>
Length Ephemeris<Frame>::AdaptiveStepParameters::length_integration_tolerance()
    const {
  return length_integration_tolerance_;
}

template<typename Frame>
Speed Ephemeris<Frame>::AdaptiveStepParameters::speed_integration_tolerance()
    const {
  return speed_integration_tolerance_;
}

template<typename Frame>
void Ephemeris<Frame>::AdaptiveStepParameters::set_max_steps(
    std::int64_t const max_steps) {
  CHECK_LT(0, max_steps);
  max_steps_ = max_steps;
}

template<typename Frame>
void Ephemeris<Frame>::AdaptiveStepParameters::set_length_integration_tolerance(
    Length const& length_integration_tolerance) {
  length_integration_tolerance_ = length_integration_tolerance;
}

template<typename Frame>
void Ephemeris<Frame>::AdaptiveStepParameters::set_speed_integration_tolerance(
    Speed const& speed_integration_tolerance) {
  speed_integration_tolerance_ = speed_integration_tolerance;
}

template<typename Frame>
void Ephemeris<Frame>::AdaptiveStepParameters::WriteToMessage(
    not_null<serialization::Ephemeris::AdaptiveStepParameters*> const message)
    const {
  integrator_->WriteToMessage(message->mutable_integrator());
  message->set_max_steps(max_steps_);
  length_integration_tolerance_.WriteToMessage(
      message->mutable_length_integration_tolerance());
  speed_integration_tolerance_.WriteToMessage(
      message->mutable_speed_integration_tolerance());
}

template<typename Frame>
typename Ephemeris<Frame>::AdaptiveStepParameters
Ephemeris<Frame>::AdaptiveStepParameters::ReadFromMessage(
    serialization::Ephemeris::AdaptiveStepParameters const& message) {
  return AdaptiveStepParameters(
      AdaptiveStepSizeIntegrator<NewtonianMotionEquation>::ReadFromMessage(
          message.integrator()),
      message.max_steps(),
      Length::ReadFromMessage(message.length_integration_tolerance()),
      Speed::ReadFromMessage(message.speed_integration_tolerance()));
}

template<typename Frame>
Ephemeris<Frame>::FixedStepParameters::FixedStepParameters(
    FixedStepSizeIntegrator<NewtonianMotionEquation> const& integrator,
    Time const& step)
    : integrator_(&integrator),
      step_(step) {
  CHECK_LT(Time(), step);
}

template<typename Frame>
inline Time const& Ephemeris<Frame>::FixedStepParameters::step() const {
  return step_;
}

template<typename Frame>
void Ephemeris<Frame>::FixedStepParameters::WriteToMessage(
    not_null<serialization::Ephemeris::FixedStepParameters*> const message)
    const {
  integrator_->WriteToMessage(message->mutable_integrator());
  step_.WriteToMessage(message->mutable_step());
}

template<typename Frame>
typename Ephemeris<Frame>::FixedStepParameters
Ephemeris<Frame>::FixedStepParameters::ReadFromMessage(
    serialization::Ephemeris::FixedStepParameters const& message) {
  return FixedStepParameters(
      FixedStepSizeIntegrator<NewtonianMotionEquation>::ReadFromMessage(
          message.integrator()),
      Time::ReadFromMessage(message.step()));
}

template<typename Frame>
Ephemeris<Frame>::Ephemeris(
    std::vector<not_null<std::unique_ptr<MassiveBody const>>>&& bodies,
    std::vector<DegreesOfFreedom<Frame>> const& initial_state,
    Instant const& initial_time,
    AccuracyParameters const& accuracy_parameters,
    FixedStepParameters const& fixed_step_parameters)
    : accuracy_parameters_(accuracy_parameters),
      fixed_step_parameters_(fixed_step_parameters) {
  CHECK(!bodies.empty());
  CHECK_EQ(bodies.size(), initial_state.size());

  IntegrationProblem<NewtonianMotionEquation> problem;
  problem.equation.compute_acceleration = [this](
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) {
    ComputeMassiveBodiesGravitationalAccelerations(t,
                                                   positions,
                                                   accelerations);
    return Status::OK;
  };

  typename NewtonianMotionEquation::SystemState& state = problem.initial_state;
  state.time = DoublePrecision<Instant>(initial_time);

  for (int i = 0; i < bodies.size(); ++i) {
    auto& body = bodies[i];
    DegreesOfFreedom<Frame> const& degrees_of_freedom = initial_state[i];

    unowned_bodies_.emplace_back(body.get());
    unowned_bodies_indices_.emplace(body.get(), i);

    auto const inserted = bodies_to_trajectories_.emplace(
                              body.get(),
                              std::make_unique<ContinuousTrajectory<Frame>>(
                                  fixed_step_parameters_.step_,
                                  accuracy_parameters_.fitting_tolerance_));
    CHECK(inserted.second);
    ContinuousTrajectory<Frame>* const trajectory =
        inserted.first->second.get();
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

  instance_ = fixed_step_parameters_.integrator_->NewInstance(
      problem,
      /*append_state=*/std::bind(
          &Ephemeris::AppendMassiveBodiesState, this, _1),
      fixed_step_parameters_.step_);
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
  absl::ReaderMutexLock l(&lock_);
  for (auto const& pair : bodies_to_trajectories_) {
    auto const& trajectory = pair.second;
    if (trajectory->empty()) {
      return true;
    }
  }
  return false;
}

template<typename Frame>
Instant Ephemeris<Frame>::t_min() const {
  absl::ReaderMutexLock l(&lock_);
  Instant t_min = bodies_to_trajectories_.begin()->second->t_min();
  for (auto const& pair : bodies_to_trajectories_) {
    auto const& trajectory = pair.second;
    t_min = std::max(t_min, trajectory->t_min());
  }
  CHECK(checkpoints_.empty() ||
        checkpoints_.front().instance->time().value >= t_min);
  return t_min;
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
  return *fixed_step_parameters_.integrator_;
}

template<typename Frame>
Status Ephemeris<Frame>::last_severe_integration_status() const {
  return last_severe_integration_status_;
}

template<typename Frame>
void Ephemeris<Frame>::ForgetBefore(Instant const& t) {
  auto it = std::upper_bound(
                checkpoints_.begin(), checkpoints_.end(), t,
                [](Instant const& left, Checkpoint const& right) {
                  // This lambda must implement a < comparison.
                  for (auto const& checkpoint : right.checkpoints) {
                    if (!checkpoint.IsAfter(left)) {
                      // The individual |checkpoint| will become invalid, so
                      // |right| <= |left|.
                      return false;
                    }
                  }
                  // All the individual checkpoints will remain valid, so
                  // |left| < |right|.
                  return true;
                });
  if (it != checkpoints_.end()) {
    CHECK_LT(t, it->instance->time().value);
  }

  for (auto& pair : bodies_to_trajectories_) {
    ContinuousTrajectory<Frame>& trajectory = *pair.second;
    trajectory.ForgetBefore(t);
  }
  checkpoints_.erase(checkpoints_.begin(), it);
}

template<typename Frame>
void Ephemeris<Frame>::Prolong(Instant const& t) {
  // Note that |t| may be before the last time that we integrated and still
  // after |t_max()|.  In this case we want to make sure that the integrator
  // makes progress.
  Instant t_final;
  Instant const instance_time = this->instance_time();
  if (t <= instance_time) {
    t_final = instance_time + fixed_step_parameters_.step_;
  } else {
    t_final = t;
  }

  // Perform the integration.  Note that we may have to iterate until |t_max()|
  // actually reaches |t| because the last series may not be fully determined
  // after the first integration.
  absl::MutexLock l(&lock_);
  while (t_max_locked() < t) {
    instance_->Solve(t_final);
    t_final += fixed_step_parameters_.step_;
  }
}

template<typename Frame>
not_null<std::unique_ptr<typename Integrator<
    typename Ephemeris<Frame>::NewtonianMotionEquation>::Instance>>
Ephemeris<Frame>::NewInstance(
    std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
    IntrinsicAccelerations const& intrinsic_accelerations,
    FixedStepParameters const& parameters) {
  IntegrationProblem<NewtonianMotionEquation> problem;

  problem.equation.compute_acceleration = [this, intrinsic_accelerations](
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) {
    if (ComputeMasslessBodiesTotalAccelerations(intrinsic_accelerations,
                                                t,
                                                positions,
                                                accelerations)) {
      return Status::OK;
    } else {
      return Status(Error::OUT_OF_RANGE, "Collision detected");
    }
  };

  CHECK(!trajectories.empty());
  Instant const trajectory_last_time = (*trajectories.begin())->last().time();
  problem.initial_state.time = DoublePrecision<Instant>(trajectory_last_time);
  for (auto const& trajectory : trajectories) {
    auto const trajectory_last = trajectory->last();
    auto const last_degrees_of_freedom = trajectory_last.degrees_of_freedom();
    CHECK_EQ(trajectory_last.time(), trajectory_last_time);
    problem.initial_state.positions.emplace_back(
        last_degrees_of_freedom.position());
    problem.initial_state.velocities.emplace_back(
        last_degrees_of_freedom.velocity());
  }

  auto const append_state =
      std::bind(&Ephemeris::AppendMasslessBodiesState, _1, trajectories);

  // The construction of the instance may evaluate the degrees of freedom of the
  // bodies.
  Prolong(trajectory_last_time + parameters.step_);

  return parameters.integrator_->NewInstance(
      problem, append_state, parameters.step_);
}

template<typename Frame>
Status Ephemeris<Frame>::FlowWithAdaptiveStep(
    not_null<DiscreteTrajectory<Frame>*> const trajectory,
    IntrinsicAcceleration intrinsic_acceleration,
    Instant const& t,
    AdaptiveStepParameters const& parameters,
    std::int64_t const max_ephemeris_steps,
    bool const last_point_only) {
  Instant const& trajectory_last_time = trajectory->last().time();
  if (trajectory_last_time == t) {
    return Status::OK;
  }

  std::vector<not_null<DiscreteTrajectory<Frame>*>> const trajectories =
      {trajectory};
  std::vector<IntrinsicAcceleration> const intrinsic_accelerations =
      {std::move(intrinsic_acceleration)};
  // The |min| is here to prevent us from spending too much time computing the
  // ephemeris.  The |max| is here to ensure that we always try to integrate
  // forward.  We use |last_state_.time.value| because this is always finite,
  // contrary to |t_max()|, which is -∞ when |empty()|.
  Instant const t_final =
      std::min(std::max(instance_time() +
                            max_ephemeris_steps * fixed_step_parameters_.step(),
                        trajectory_last_time + fixed_step_parameters_.step()),
               t);
  Prolong(t_final);

  IntegrationProblem<NewtonianMotionEquation> problem;
  problem.equation.compute_acceleration = [this, &intrinsic_accelerations](
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) {
    if (ComputeMasslessBodiesTotalAccelerations(intrinsic_accelerations,
                                                t,
                                                positions,
                                                accelerations)) {
      return Status::OK;
    } else {
      return Status(Error::OUT_OF_RANGE, "Collision detected");
    }
  };

  auto const trajectory_last = trajectory->last();
  auto const last_degrees_of_freedom = trajectory_last.degrees_of_freedom();
  problem.initial_state = {{last_degrees_of_freedom.position()},
                           {last_degrees_of_freedom.velocity()},
                           trajectory_last.time()};

  typename AdaptiveStepSizeIntegrator<NewtonianMotionEquation>::Parameters const
      integrator_parameters(
          /*first_time_step=*/t_final - problem.initial_state.time.value,
          /*safety_factor=*/0.9,
          parameters.max_steps_,
          /*last_step_is_exact=*/true);
  CHECK_GT(integrator_parameters.first_time_step, 0 * Second)
      << "Flow back to the future: " << t_final
      << " <= " << problem.initial_state.time.value;
  auto const tolerance_to_error_ratio =
      std::bind(&Ephemeris<Frame>::ToleranceToErrorRatio,
                std::cref(parameters.length_integration_tolerance_),
                std::cref(parameters.speed_integration_tolerance_),
                _1, _2);

  typename AdaptiveStepSizeIntegrator<NewtonianMotionEquation>::AppendState
      append_state;
  typename NewtonianMotionEquation::SystemState last_state;
  if (last_point_only) {
    append_state = [&last_state](
        typename NewtonianMotionEquation::SystemState const& state) {
      last_state = state;
    };
  } else {
    append_state = std::bind(
        &Ephemeris::AppendMasslessBodiesState, _1, std::cref(trajectories));
  }

  auto const instance =
      parameters.integrator_->NewInstance(problem,
                                          append_state,
                                          tolerance_to_error_ratio,
                                          integrator_parameters);
  auto status = instance->Solve(t_final);

  // We probably don't care if the vessel gets too close to the singularity, as
  // we only use this integrator for the future.  So we swallow the error.
  // TODO(phl): Is this the right thing to do long term?
  if (status.error() == Error::OUT_OF_RANGE) {
    status = Status::OK;
  }

  if (last_point_only) {
    AppendMasslessBodiesState(last_state, trajectories);
  }

  // TODO(egg): when we have events in trajectories, we should add a singularity
  // event at the end if the outcome indicates a singularity
  // (|VanishingStepSize|).  We should not have an event on the trajectory if
  // |ReachedMaximalStepCount|, since that is not a physical property, but
  // rather a self-imposed constraint.
  if (!status.ok() || t_final == t) {
    return status;
  } else {
    return Status(Error::DEADLINE_EXCEEDED,
                  "Couldn't reach " + DebugString(t_final) + ", stopping at " +
                      DebugString(t));
  }
}

template<typename Frame>
Status Ephemeris<Frame>::FlowWithAdaptiveStepGeneralized(
      not_null<DiscreteTrajectory<Frame>*> trajectory,
      GeneralIntrinsicAcceleration intrinsic_acceleration,
      Instant const& t,
      AdaptiveStepParameters const& parameters,
      std::int64_t max_ephemeris_steps,
      bool last_point_only) {
  // TODO(phl): FlowWithAdaptiveStepGeneralized and FlowWithAdaptiveStep are
  // very similar.  Factor that which can be factored.
  using ODE = ExplicitSecondOrderOrdinaryDifferentialEquation<Position<Frame>>;
  Instant const& trajectory_last_time = trajectory->last().time();
  if (trajectory_last_time == t) {
    return Status::OK;
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
  Prolong(t_final);

  IntegrationProblem<ODE> problem;
  problem.equation.compute_acceleration =
      [this, &intrinsic_acceleration](
          Instant const& t,
          std::vector<Position<Frame>> const& positions,
          std::vector<Velocity<Frame>> const& velocities,
          std::vector<Vector<Acceleration, Frame>>& accelerations) {
        if (!ComputeMasslessBodiesGravitationalAccelerations(
                t, positions, accelerations)) {
          return Status(Error::OUT_OF_RANGE, "Collision detected");
        }
        accelerations[0] +=
            intrinsic_acceleration(t, {positions[0], velocities[0]});
        return Status::OK;
      };

  auto const trajectory_last = trajectory->last();
  auto const last_degrees_of_freedom = trajectory_last.degrees_of_freedom();
  problem.initial_state = {{last_degrees_of_freedom.position()},
                           {last_degrees_of_freedom.velocity()},
                           trajectory_last.time()};

  typename AdaptiveStepSizeIntegrator<ODE>::Parameters const
      integrator_parameters(
          /*first_time_step=*/t_final - problem.initial_state.time.value,
          /*safety_factor=*/0.9,
          parameters.max_steps_,
          /*last_step_is_exact=*/true);
  CHECK_GT(integrator_parameters.first_time_step, 0 * Second)
      << "Flow back to the future: " << t_final
      << " <= " << problem.initial_state.time.value;
  auto const tolerance_to_error_ratio =
      std::bind(&Ephemeris<Frame>::ToleranceToErrorRatio,
                std::cref(parameters.length_integration_tolerance_),
                std::cref(parameters.speed_integration_tolerance_),
                _1, _2);

  typename AdaptiveStepSizeIntegrator<ODE>::AppendState append_state;
  typename ODE::SystemState last_state;
  if (last_point_only) {
    append_state = [&last_state](typename ODE::SystemState const& state) {
      last_state = state;
    };
  } else {
    append_state = std::bind(
        &Ephemeris::AppendMasslessBodiesState, _1, std::cref(trajectories));
  }

  auto const instance =
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
          Fine1987RKNG34, Position<Frame>>().NewInstance(
              problem,
              append_state,
              tolerance_to_error_ratio,
              integrator_parameters);
  auto status = instance->Solve(t_final);

  // We probably don't care if the vessel gets too close to the singularity, as
  // we only use this integrator for the future.  So we swallow the error.
  // TODO(phl): Is this the right thing to do long term?
  if (status.error() == Error::OUT_OF_RANGE) {
    status = Status::OK;
  }

  if (last_point_only) {
    AppendMasslessBodiesState(last_state, trajectories);
  }

  // TODO(egg): when we have events in trajectories, we should add a singularity
  // event at the end if the outcome indicates a singularity
  // (|VanishingStepSize|).  We should not have an event on the trajectory if
  // |ReachedMaximalStepCount|, since that is not a physical property, but
  // rather a self-imposed constraint.
  if (!status.ok() || t_final == t) {
    return status;
  } else {
    return Status(Error::DEADLINE_EXCEEDED,
                  "Couldn't reach " + DebugString(t_final) + ", stopping at " +
                      DebugString(t));
  }
}

template<typename Frame>
Status Ephemeris<Frame>::FlowWithFixedStep(
    Instant const& t,
    typename Integrator<NewtonianMotionEquation>::Instance& instance) {
  if (empty() || t > t_max()) {
    Prolong(t);
  }

  return instance.Solve(t);
}

template<typename Frame>
Vector<Acceleration, Frame> Ephemeris<Frame>::
ComputeGravitationalAccelerationOnMasslessBody(
    Position<Frame> const& position,
    Instant const& t) const {
  std::vector<Vector<Acceleration, Frame>> accelerations(1);
  ComputeMasslessBodiesGravitationalAccelerations(t, {position}, accelerations);

  return accelerations[0];
}

template<typename Frame>
Vector<Acceleration, Frame> Ephemeris<Frame>::
ComputeGravitationalAccelerationOnMasslessBody(
    not_null<DiscreteTrajectory<Frame>*> const trajectory,
    Instant const& t) const {
  auto const it = trajectory->Find(t);
  DegreesOfFreedom<Frame> const& degrees_of_freedom = it.degrees_of_freedom();
  return ComputeGravitationalAccelerationOnMasslessBody(
             degrees_of_freedom.position(), t);
}

template<typename Frame>
Vector<Acceleration, Frame> Ephemeris<Frame>::
ComputeGravitationalAccelerationOnMassiveBody(
    not_null<MassiveBody const*> const body,
    Instant const& t) const {
  bool const body_is_oblate = body->is_oblate();

  std::vector<Position<Frame>> positions;
  std::vector<Vector<Acceleration, Frame>> accelerations(bodies_.size());
  int b1 = -1;

  // Evaluate the |positions|.
  positions.reserve(bodies_.size());
  for (int b = 0; b < bodies_.size(); ++b) {
    auto const& current_body = bodies_[b];
    auto const& current_body_trajectory = trajectories_[b];
    if (current_body.get() == body) {
      CHECK_EQ(-1, b1);
      b1 = b;
    }
    positions.push_back(current_body_trajectory->EvaluatePosition(t));
  }
  CHECK_LE(0, b1);

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
void Ephemeris<Frame>::ComputeApsides(not_null<MassiveBody const*> const body1,
                                      not_null<MassiveBody const*> const body2,
                                      DiscreteTrajectory<Frame>& apoapsides1,
                                      DiscreteTrajectory<Frame>& periapsides1,
                                      DiscreteTrajectory<Frame>& apoapsides2,
                                      DiscreteTrajectory<Frame>& periapsides2) {
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
        body1_trajectory->EvaluateDegreesOfFreedom(t);
    DegreesOfFreedom<Frame> const body2_degrees_of_freedom =
        body2_trajectory->EvaluateDegreesOfFreedom(t);
    RelativeDegreesOfFreedom<Frame> const relative =
        body1_degrees_of_freedom - body2_degrees_of_freedom;
    return 2.0 * InnerProduct(relative.displacement(), relative.velocity());
  };

  std::optional<Instant> previous_time;
  std::optional<Variation<Square<Length>>> previous_squared_distance_derivative;

  for (Instant time = t_min();
       time <= t_max();
       time += fixed_step_parameters_.step()) {
    Variation<Square<Length>> const squared_distance_derivative =
        evaluate_square_distance_derivative(time);
    if (previous_squared_distance_derivative &&
        Sign(squared_distance_derivative) !=
            Sign(*previous_squared_distance_derivative)) {
      CHECK(previous_time);

      // The derivative of |squared_distance| changed sign.  Find its zero by
      // bisection, this is the time of the apsis.  Then compute the apsis and
      // append it to one of the output trajectories.
      Instant const apsis_time = Bisect(evaluate_square_distance_derivative,
                                        *previous_time,
                                        time);
      DegreesOfFreedom<Frame> const apsis1_degrees_of_freedom =
          body1_trajectory->EvaluateDegreesOfFreedom(apsis_time);
      DegreesOfFreedom<Frame> const apsis2_degrees_of_freedom =
          body2_trajectory->EvaluateDegreesOfFreedom(apsis_time);
      if (Sign(squared_distance_derivative).Negative()) {
        apoapsides1.Append(apsis_time, apsis1_degrees_of_freedom);
        apoapsides2.Append(apsis_time, apsis2_degrees_of_freedom);
      } else {
        periapsides1.Append(apsis_time, apsis1_degrees_of_freedom);
        periapsides2.Append(apsis_time, apsis2_degrees_of_freedom);
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
  // The bodies are serialized in the order in which they were given at
  // construction.
  for (auto const& unowned_body : unowned_bodies_) {
    unowned_body->WriteToMessage(message->add_body());
  }
  // The trajectories are serialized in the order resulting from the separation
  // between oblate and spherical bodies.
  if (checkpoints_.empty()) {
    for (auto const& trajectory : trajectories_) {
      trajectory->WriteToMessage(message->add_trajectory());
    }
    instance_->WriteToMessage(message->mutable_instance());
  } else {
    auto const& checkpoints = checkpoints_.front().checkpoints;
    CHECK_EQ(trajectories_.size(), checkpoints.size());
    for (int i = 0; i < trajectories_.size(); ++i) {
      trajectories_[i]->WriteToMessage(message->add_trajectory(),
                                       checkpoints[i]);
    }
    checkpoints_.front().instance->WriteToMessage(
        message->mutable_instance());
    t_max().WriteToMessage(message->mutable_t_max());
  }
  fixed_step_parameters_.WriteToMessage(
      message->mutable_fixed_step_parameters());
  accuracy_parameters_.WriteToMessage(
      message->mutable_accuracy_parameters());
  LOG(INFO) << NAMED(message->SpaceUsed());
  LOG(INFO) << NAMED(message->ByteSize());
}

template<typename Frame>
not_null<std::unique_ptr<Ephemeris<Frame>>> Ephemeris<Frame>::ReadFromMessage(
    serialization::Ephemeris const& message) {
  bool const is_pre_ἐρατοσθένης = !message.has_accuracy_parameters();
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  for (auto const& body : message.body()) {
    bodies.push_back(MassiveBody::ReadFromMessage(body));
  }

  AccuracyParameters accuracy_parameters(
      pre_ἐρατοσθένης_default_ephemeris_fitting_tolerance);
  if (!is_pre_ἐρατοσθένης) {
    accuracy_parameters =
        AccuracyParameters::ReadFromMessage(message.accuracy_parameters());
  }
  FixedStepParameters const fixed_step_parameters =
      FixedStepParameters::ReadFromMessage(message.fixed_step_parameters());

  // Dummy initial state and time.  We'll overwrite them later.
  std::vector<DegreesOfFreedom<Frame>> const initial_state(
      bodies.size(),
      DegreesOfFreedom<Frame>(Position<Frame>(), Velocity<Frame>()));
  Instant const initial_time;
  auto ephemeris = make_not_null_unique<Ephemeris<Frame>>(
                       std::move(bodies),
                       initial_state,
                       initial_time,
                       accuracy_parameters,
                       fixed_step_parameters);

  NewtonianMotionEquation equation;
  equation.compute_acceleration = [ephemeris = ephemeris.get()](
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) {
    ephemeris->ComputeMassiveBodiesGravitationalAccelerations(t,
                                                              positions,
                                                              accelerations);
    return Status::OK;
  };

  ephemeris->instance_ =
      FixedStepSizeIntegrator<NewtonianMotionEquation>::Instance::
      ReadFromMessage(
          message.instance(),
          equation,
          /*append_state=*/std::bind(
              &Ephemeris::AppendMassiveBodiesState, ephemeris.get(), _1));

  int index = 0;
  ephemeris->bodies_to_trajectories_.clear();
  ephemeris->trajectories_.clear();
  for (auto const& trajectory : message.trajectory()) {
    not_null<MassiveBody const*> const body = ephemeris->bodies_[index].get();
    not_null<std::unique_ptr<ContinuousTrajectory<Frame>>>
        deserialized_trajectory =
            ContinuousTrajectory<Frame>::ReadFromMessage(trajectory);
    ephemeris->trajectories_.push_back(deserialized_trajectory.get());
    ephemeris->bodies_to_trajectories_.emplace(
        body, std::move(deserialized_trajectory));
    ++index;
  }
  if (message.has_t_max()) {
    ephemeris->checkpoints_.push_back(ephemeris->GetCheckpoint());
    ephemeris->Prolong(Instant::ReadFromMessage(message.t_max()));
  }
  return ephemeris;
}

template<typename Frame>
Ephemeris<Frame>::Ephemeris(
    FixedStepSizeIntegrator<
        typename Ephemeris<Frame>::NewtonianMotionEquation> const& integrator)
    : accuracy_parameters_(pre_ἐρατοσθένης_default_ephemeris_fitting_tolerance),
      fixed_step_parameters_(integrator, 1 * Second) {}

template<typename Frame>
void Ephemeris<Frame>::AppendMassiveBodiesState(
    typename NewtonianMotionEquation::SystemState const& state) {
  int index = 0;
  for (int i = 0; i < trajectories_.size(); ++i) {
    auto const& trajectory = trajectories_[i];
    auto const status = trajectory->Append(
        state.time.value,
        DegreesOfFreedom<Frame>(state.positions[index].value,
                                state.velocities[index].value));

    // Handle the apocalypse.
    if (!status.ok()) {
      last_severe_integration_status_ =
          Status(status.error(),
                 "Error extending trajectory for " + bodies_[i]->name() + ". " +
                     status.message());
      LOG(ERROR) << "New Apocalypse: " << last_severe_integration_status_;
    }

    ++index;
  }

  // Record an intermediate state if we haven't done so for too long.
  CHECK(!trajectories_.empty());
  Instant const t_last_intermediate_state =
      checkpoints_.empty()
          ? astronomy::InfinitePast
          : checkpoints_.back().instance->time().value;
  if (t_max_locked() - t_last_intermediate_state >
      max_time_between_checkpoints) {
    checkpoints_.push_back(GetCheckpoint());
  }
}

template<typename Frame>
void Ephemeris<Frame>::AppendMasslessBodiesState(
    typename NewtonianMotionEquation::SystemState const& state,
    std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories) {
  int index = 0;
  for (auto& trajectory : trajectories) {
    trajectory->Append(
        state.time.value,
        DegreesOfFreedom<Frame>(state.positions[index].value,
                                state.velocities[index].value));
    ++index;
  }
}

template<typename Frame>
typename Ephemeris<Frame>::Checkpoint Ephemeris<Frame>::GetCheckpoint() {
  std::vector<typename ContinuousTrajectory<Frame>::Checkpoint> checkpoints;
  for (auto const& trajectory : trajectories_) {
    checkpoints.push_back(trajectory->GetCheckpoint());
  }
  return Checkpoint({instance_->Clone(), checkpoints});
}

template<typename Frame>
Instant Ephemeris<Frame>::t_max_locked() const {
  Instant t_max = bodies_to_trajectories_.begin()->second->t_max();
  for (auto const& pair : bodies_to_trajectories_) {
    auto const& trajectory = pair.second;
    t_max = std::min(t_max, trajectory->t_max());
  }
  // Here we may have a checkpoint after |t_max| if the checkpointed state was
  // not yet incorporated in a series.
  return t_max;
}

template<typename Frame>
Instant Ephemeris<Frame>::instance_time() const {
  absl::ReaderMutexLock l(&lock_);
  return instance_->time().value;
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

    // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
    // sive corporum duorum actiones in se mutuo semper esse æquales &
    // in partes contrarias dirigi.
    auto const μ2_over_Δq³ = μ2 * one_over_Δq³;
    acceleration_on_b1 -= Δq * μ2_over_Δq³;

    if (body1_is_oblate || body2_is_oblate) {
      if (body1_is_oblate) {
        Vector<Quotient<Acceleration,
                        GravitationalParameter>, Frame> const
            degree_2_zonal_effect1 =
                geopotentials[b1].GeneralSphericalHarmonicsAcceleration(
                    t,
                    -Δq,
                    Δq_norm,
                    Δq²,
                    one_over_Δq³);
        acceleration_on_b1 -= μ2 * degree_2_zonal_effect1;
        acceleration_on_b2 += μ1 * degree_2_zonal_effect1;
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
bool Ephemeris<Frame>::
ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies(
    Instant const& t,
    MassiveBody const& body1,
    std::size_t const b1,
    std::vector<Position<Frame>> const& positions,
    std::vector<Vector<Acceleration, Frame>>& accelerations) const {
  GravitationalParameter const& μ1 = body1.gravitational_parameter();
  Position<Frame> const position1 = trajectories_[b1]->EvaluatePosition(t);
  Length const body1_mean_radius = body1.mean_radius();
  bool ok = true;

  for (std::size_t b2 = 0; b2 < positions.size(); ++b2) {
    // A vector from the center of |b2| to the center of |b1|.
    Displacement<Frame> const Δq = position1 - positions[b2];

    Square<Length> const Δq² = Δq.Norm²();
    Length const Δq_norm = Sqrt(Δq²);
    ok &= Δq_norm > body1_mean_radius;

    Exponentiation<Length, -3> const one_over_Δq³ = Δq_norm / (Δq² * Δq²);

    auto const μ1_over_Δq³ = μ1 * one_over_Δq³;
    accelerations[b2] += Δq * μ1_over_Δq³;

    if (body1_is_oblate) {
      Vector<Quotient<Acceleration,
                      GravitationalParameter>, Frame> const
          degree_2_zonal_effect1 =
              geopotentials_[b1].GeneralSphericalHarmonicsAcceleration(
                  t,
                  -Δq,
                  Δq_norm,
                  Δq²,
                  one_over_Δq³);
      accelerations[b2] += μ1 * degree_2_zonal_effect1;
    }
  }
  return ok;
}

template<typename Frame>
void Ephemeris<Frame>::ComputeMassiveBodiesGravitationalAccelerations(
    Instant const& t,
    std::vector<Position<Frame>> const& positions,
    std::vector<Vector<Acceleration, Frame>>& accelerations) const {
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
}

template<typename Frame>
bool Ephemeris<Frame>::ComputeMasslessBodiesGravitationalAccelerations(
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) const {
  CHECK_EQ(positions.size(), accelerations.size());
  accelerations.assign(accelerations.size(), Vector<Acceleration, Frame>());
  bool ok = true;

  absl::ReaderMutexLock l(&lock_);
  for (std::size_t b1 = 0; b1 < number_of_oblate_bodies_; ++b1) {
    MassiveBody const& body1 = *bodies_[b1];
    ok &= ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies<
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
    ok &= ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies<
        /*body1_is_oblate=*/false>(
        t,
        body1, b1,
        positions,
        accelerations);
  }
  return ok;
}

template<typename Frame>
bool Ephemeris<Frame>::ComputeMasslessBodiesTotalAccelerations(
    IntrinsicAccelerations const& intrinsic_accelerations,
    Instant const& t,
    std::vector<Position<Frame>> const& positions,
    std::vector<Vector<Acceleration, Frame>>& accelerations) const {
  // First, the acceleration due to the gravitational field of the
  // massive bodies.
  bool const ok = ComputeMasslessBodiesGravitationalAccelerations(
      t, positions, accelerations);

  // Then, the intrinsic accelerations, if any.
  if (!intrinsic_accelerations.empty()) {
    for (int i = 0; i < intrinsic_accelerations.size(); ++i) {
      auto const intrinsic_acceleration = intrinsic_accelerations[i];
      if (intrinsic_acceleration != nullptr) {
        accelerations[i] += intrinsic_acceleration(t);
      }
    }
  }
  return ok;
}

template<typename Frame>
double Ephemeris<Frame>::ToleranceToErrorRatio(
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance,
    Time const& current_step_size,
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
