
#ifndef PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
#define PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_

#include <experimental/optional>
#include <functional>
#include <limits>
#include <vector>

#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/double_precision.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace internal_ordinary_differential_equations {

using base::Error;
using base::not_null;
using base::Status;
using geometry::Instant;
using numerics::DoublePrecision;
using quantities::Time;
using quantities::Variation;

// A differential equation of the form q″ = f(q, t).
// |Position| is the type of q.
template<typename Position_>
struct SpecialSecondOrderDifferentialEquation {
  using Position = Position_;
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
               std::vector<Acceleration>& accelerations)>;

  struct SystemState {
    std::vector<DoublePrecision<Position>> positions;
    std::vector<DoublePrecision<Velocity>> velocities;
    DoublePrecision<Instant> time;

    void WriteToMessage(
        not_null<serialization::SystemState*> const message) const;
    static SystemState ReadFromMessage(
        serialization::SystemState const& message);
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

// An initial value problem.
template<typename ODE>
struct IntegrationProblem {
  ODE equation;
  typename ODE::SystemState const* initial_state;
};

// An opaque object for holding the state during the integration of a problem.
struct IntegrationInstance {
  template<typename ODE>
  using AppendState =
      std::function<void(typename ODE::SystemState const& state)>;
  virtual ~IntegrationInstance() = default;  // Makes the type polymorphic.
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
  // This number must be in ]0, 1[.  Higher values increase the chance of step
  // rejection, lower values yield smaller steps.
  double safety_factor;
  // This functor is called at each step, with the |current_step_size| used by
  // the integrator and the estimated |error| on that step.  It returns the
  // ratio of a tolerance to some norm of the error.  The step is recomputed
  // with a smaller step size if the result is less than 1, and accepted
  // otherwise.
  // In both cases, the new step size is chosen so as to try and make the result
  // of the next call to |tolerance_to_error_ratio| close to |safety_factor|.
  ToleranceToErrorRatio tolerance_to_error_ratio;
  // Integration will stop after |*max_steps| even if it has not reached
  // |t_final|.
  std::int64_t max_steps = std::numeric_limits<std::int64_t>::max();
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
  using ODE = DifferentialEquation;
  // The last call to |append_state| has a |state.time.value| equal to the
  // unique |Instant| of the form |problem.t_final + n * step| in
  // ]t_final - step, t_final].  |append_state| will be called with
  // |state.time.values|s at intervals differing from |step| by at most one ULP.
  virtual void Solve(Instant const& t_final,
                     IntegrationInstance& instance) const = 0;

  virtual not_null<std::unique_ptr<IntegrationInstance>> NewInstance(
    IntegrationProblem<ODE> const& problem,
    IntegrationInstance::AppendState<ODE> append_state,
    Time const& step) const = 0;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> const message) const;
  static FixedStepSizeIntegrator const& ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message);

 protected:
  explicit FixedStepSizeIntegrator(
      serialization::FixedStepSizeIntegrator::Kind const kind);

 private:
  serialization::FixedStepSizeIntegrator::Kind const kind_;
};

}  // namespace internal_ordinary_differential_equations

// The |Solve| function below exclusively returns one of the following statuses.
namespace termination_condition {
constexpr base::Error Done = base::Error::OK;
// The integration may be retried with the same arguments and progress will
// happen.
constexpr base::Error ReachedMaximalStepCount = base::Error::ABORTED;
// A singularity.
constexpr base::Error VanishingStepSize = base::Error::FAILED_PRECONDITION;
}  // namespace termination_condition

namespace internal_ordinary_differential_equations {

// An integrator using a fixed step size.
template<typename DifferentialEquation>
class AdaptiveStepSizeIntegrator : public Integrator<DifferentialEquation> {
 public:
  using ODE = DifferentialEquation;
  // The last call to |append_state| will have |state.time.value == t_final|.
  virtual Status Solve(Instant const& t_final,
                       IntegrationInstance& instance) const = 0;

  virtual not_null<std::unique_ptr<IntegrationInstance>> NewInstance(
    IntegrationProblem<ODE> const& problem,
    IntegrationInstance::AppendState<ODE> append_state,
    AdaptiveStepSize<ODE> const& adaptive_step_size) const = 0;

  void WriteToMessage(
      not_null<serialization::AdaptiveStepSizeIntegrator*> const message) const;
  static AdaptiveStepSizeIntegrator const& ReadFromMessage(
      serialization::AdaptiveStepSizeIntegrator const& message);

 protected:
  explicit AdaptiveStepSizeIntegrator(
      serialization::AdaptiveStepSizeIntegrator::Kind const kind);

 private:
  serialization::AdaptiveStepSizeIntegrator::Kind const kind_;
};

}  // namespace internal_ordinary_differential_equations

using internal_ordinary_differential_equations::AdaptiveStepSizeIntegrator;
using internal_ordinary_differential_equations::FixedStepSizeIntegrator;
using internal_ordinary_differential_equations::IntegrationInstance;
using internal_ordinary_differential_equations::IntegrationProblem;
using internal_ordinary_differential_equations::Integrator;
using internal_ordinary_differential_equations::SpecialSecondOrderDifferentialEquation;

}  // namespace integrators
}  // namespace principia

#include "integrators/ordinary_differential_equations_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
