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

// A differential equation of the form q″ = f(q, t).
// |Position| is the type of q.
template<typename Position>
struct SpecialSecondOrderDifferentialEquation {
  // The type of Δq.
  using Displacement = Difference<Position>;
  // The type of q′.
  using Velocity = Variation<Position>;
  // The type of q″.
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
  // A functor that computes f(q, t) and stores it in |*accelerations|.
  // This functor must be called with |accelerations->size()| equal to
  // |positions->size()|, but there is no requirement on the values in
  // |*acceleration|.
  RightHandSideComputation compute_acceleration;
};

// An initial value problem, together with a final time for the solution
// and a callback for processing solution points.
template<typename ODE>
struct IntegrationProblem {
  ODE equation;
  typename ODE::SystemState const* initial_state;
  Instant t_final;
  std::function<void(typename ODE::SystemState const& state)> append_state;
};

// Settings for for adaptive step size integration.
template<typename ODE>
struct AdaptiveStepSize {
  using ToleranceToErrorRatio =
      std::function<
          double(Time const& current_step_size,
                 typename ODE::SystemStateError const& error)>;
  // The first time step tried by the integrator. It must have the same sign as
  // |problem.t_final - initial_state.time.value|.
  Time first_time_step;
  // This number must be in ]0, 1[.
  double safety_factor;
  // This functor is called at each step, with the |current_step_size| used by
  // the integrator and the estimated |error| on that step.  It returns the
  // ratio of a tolerance to some norm of the error.  The step is recomputed
  // with a smaller step size if the result is less than 1, and accepted
  // otherwise.
  // In both cases, the new step size is chosen so as to try and make the result
  // of the next call to |tolerance_to_error_ratio| close to |safety_factor|.
  ToleranceToErrorRatio tolerance_to_error_ratio;
};

// A base class for integrators.
template<typename DifferentialEquation>
class Integrator {
 public:
  using ODE = DifferentialEquation;
};

// An integrator using a fixed step size.
template<typename DifferentialEquation>
class FixedStepSizeIntegrator : public Integrator<DifferentialEquation> {
 public:
  // The last call to |problem.append_state| has a |state.time.value| equal to
  // the unique |Instant| of the form |problem.t_final + n * step| in
  // [problem.t_final, problem.t_final + step[.
  // |problem.append_state| will be called with |state.time.values|s at
  // intervals differing from |step| by at most one ULP.
  virtual void Solve(IntegrationProblem<ODE> const& problem,
                     Time const& step) const = 0;
};

// An integrator using a fixed step size.
template<typename DifferentialEquation>
class AdaptiveStepSizeIntegrator : public Integrator<DifferentialEquation> {
 public:
  // The last call to |problem.append_state| will have
  // |state.time.value == problem.t_final|.
  virtual void Solve(IntegrationProblem<ODE> const& problem,
                     AdaptiveStepSize<ODE> const& adaptive_step_size) const = 0;
};

}  // namespace integrators
}  // namespace principia
