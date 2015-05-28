#pragma once

#include <functional>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/double_precision.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using base::not_null;
using geometry::Instant;
using numerics::DoublePrecision;
using quantities::Time;
using quantities::Variation;

namespace integrators {

template<typename Position>
struct SpecialSecondOrderDifferentialEquation {
  using Displacement = Difference<Position>;
  using Velocity = Variation<Position>;
  using Acceleration = Variation<Velocity>;
  using RightHandSideComputation =
      std::function<
          void(Instant const& t,
               std::vector<Position> const& positions,
               not_null<std::vector<Acceleration>*> const accelerations)>;
  struct SystemState {
    std::vector<DoublePrecision<Position>> positions;
    std::vector<DoublePrecision<Velocity>> velocities;
    DoublePrecision<Instant> time;
  };
  struct SystemStateError {
    std::vector<Displacement> position_error;
    std::vector<Velocity> velocity_error;
  };
  RightHandSideComputation compute_acceleration;
};

template<typename ODE>
struct IntegrationProblem {
  ODE equation;
  typename ODE::SystemState const* initial_state;
  Instant t_final;
  std::function<void(typename ODE::SystemState const& state)> append_state;
};

template<typename ODE>
struct AdaptiveStepSize {
  using ToleranceToErrorRatio =
      std::function<
          double(Time const& current_step_size,
                 typename ODE::SystemStateError const& error)>;
  Time first_time_step;
  double safety_factor;
  ToleranceToErrorRatio tolerance_to_error_ratio;
};

template<typename DifferentialEquation>
class Integrator {
 public:
  using ODE = DifferentialEquation;
};

template<typename DifferentialEquation>
class FixedStepSizeIntegrator : public Integrator<DifferentialEquation> {
 public:
  virtual void Solve(IntegrationProblem<ODE> const& problem,
                     Time const& step) const = 0;
};

template<typename DifferentialEquation>
class AdaptiveStepSizeIntegrator : public Integrator<DifferentialEquation> {
 public:
  virtual void Solve(IntegrationProblem<ODE> const& problem,
                     AdaptiveStepSize<ODE> const& adaptive_step_size) const = 0;
};

}  // namespace integrators
}  // namespace principia
