#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#define PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_

#include <functional>
#include <string>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace _integrators {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_instant;
using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_quantities;

// A base class for integrators.
template<typename ODE_>
class Integrator {
 public:
  using ODE = ODE_;
  using IndependentVariable = typename ODE::IndependentVariable;
  using AppendState =
      std::function<void(typename ODE::State const& state)>;

  // An object for holding the integrator state during the integration of a
  // problem.
  class Instance {
   public:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state);
    virtual ~Instance() = default;

    // The subclass must document the time (or independent variable) passed to
    // the last call to |append_state|.
    virtual absl::Status Solve(IndependentVariable const& s_final) = 0;

    // The equation integrated by this instance.
    ODE const& equation() const;

    // The |AppendState| function.
    AppendState const& append_state() const;

    // The last instant integrated by this instance.
    DoublePrecision<Instant> const& time() const;

    // The last state integrated by this instance.
    typename ODE::State const& state() const;
    typename ODE::State& state();

    // Performs a copy of this object.
    virtual not_null<std::unique_ptr<Instance>> Clone() const = 0;

    // |ReadFromMessage| is specific to each subclass because of the functions.
    virtual void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const;

   protected:
    // For testing.
    Instance();

    // We make the data members protected because they need to be easily
    // accessible by subclasses.
    ODE const equation_;
    typename ODE::State current_state_;
    AppendState const append_state_;
  };

  virtual ~Integrator() = default;
};

// An integrator using a fixed step size.
template<typename ODE_>
class FixedStepSizeIntegrator : public Integrator<ODE_> {
 public:
  using ODE = ODE_;
  using typename Integrator<ODE>::AppendState;

  // The last call to |append_state| has a |state.time.value| equal to the
  // unique |Instant| of the form |t_final + n * step| in
  // ]t_final - step, t_final].  |append_state| will be called with
  // |state.time.values|s at intervals differing from |step| by at most one ULP.
  class Instance : public Integrator<ODE>::Instance {
   public:
    // The time step used by this instance.
    Time const& step() const;

    // The integrator corresponding to this instance.
    virtual FixedStepSizeIntegrator const& integrator() const = 0;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
    template<typename S = typename ODE::State,
             typename = std::enable_if_t<is_serializable_v<S>>>
    static not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
    ReadFromMessage(serialization::IntegratorInstance const& message,
                    ODE const& equation,
                    AppendState const& append_state);

   protected:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state,
             Time const& step);

    Time const step_;
  };

  // The factory function for |Instance|, above.  It ensures that the instance
  // has a back-pointer to its integrator.
  virtual not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
  NewInstance(InitialValueProblem<ODE> const& problem,
              AppendState const& append_state,
              Time const& step) const = 0;

  virtual void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const = 0;
  static FixedStepSizeIntegrator const& ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message);

 protected:
  FixedStepSizeIntegrator() = default;
};

template<typename Equation>
FixedStepSizeIntegrator<Equation> const&
ParseFixedStepSizeIntegrator(std::string const& integrator_kind);

// An integrator using an adaptive step size.
template<typename ODE_>
class AdaptiveStepSizeIntegrator : public Integrator<ODE_> {
 public:
  using ODE = ODE_;
  using IndependentVariableDifference =
      typename ODE::IndependentVariableDifference;
  using typename Integrator<ODE>::AppendState;

  // This functor is called at each step, with the |current_step_size| used by
  // the integrator and the estimated |error| on that step.  It returns the
  // ratio of a tolerance to some norm of the error.  The step is recomputed
  // with a smaller step size if the result is less than 1, and accepted
  // otherwise.
  // In both cases, the new step size is chosen so as to try and make the
  // result of the next call to this functor close to |safety_factor|.
  using ToleranceToErrorRatio = std::function<double(
      IndependentVariableDifference const& current_step_size,
      typename ODE::State const& state,
      typename ODE::State::Error const& error)>;

  struct Parameters final {
    Parameters(IndependentVariableDifference const& first_step,
               double safety_factor,
               std::int64_t max_steps,
               bool last_step_is_exact);

    // |max_steps| is infinite and the last step is exact.
    Parameters(IndependentVariableDifference const& first_step,
               double safety_factor);

    void WriteToMessage(
        not_null<serialization::AdaptiveStepSizeIntegratorInstance::
                     Parameters*> message) const;
    static Parameters ReadFromMessage(
        serialization::AdaptiveStepSizeIntegratorInstance::Parameters const&
            message);

    // The first step tried by the integrator. It must have the same sign as
    // |problem.t_final - initial_state.time.value|.
    IndependentVariableDifference const first_step;
    // This number must be in ]0, 1[.  Higher values increase the chance of step
    // rejection, lower values yield smaller steps.
    double const safety_factor;
    // Integration will stop after |*max_steps| even if it has not reached
    // |t_final|.
    std::int64_t const max_steps;
    // If true, the he last call to |append_state| has
    // |state.time.value == t_final| (unless |max_steps| is reached).  Otherwise
    // it may have |state.time.value < t_final|.
    bool const last_step_is_exact;
  };

  // The last call to |append_state| will have |state.time.value == t_final|.
  class Instance : public Integrator<ODE>::Instance {
   public:
    // The integrator corresponding to this instance.
    virtual AdaptiveStepSizeIntegrator const& integrator() const = 0;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
    template<typename S = typename ODE::State,
             typename = std::enable_if_t<is_serializable_v<S>>>
    static not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
    ReadFromMessage(serialization::IntegratorInstance const& message,
                    ODE const& equation,
                    AppendState const& append_state,
                    ToleranceToErrorRatio const& tolerance_to_error_ratio);

   protected:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state,
             ToleranceToErrorRatio tolerance_to_error_ratio,
             Parameters const& parameters,
             IndependentVariableDifference const& step,
             bool first_use);

    ToleranceToErrorRatio const tolerance_to_error_ratio_;
    Parameters const parameters_;
    IndependentVariableDifference step_;
    bool first_use_;
  };

  // The factory function for |Instance|, above.  It ensures that the instance
  // has a back-pointer to its integrator.
  virtual not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
  NewInstance(InitialValueProblem<ODE> const& problem,
              typename Integrator<ODE>::AppendState const& append_state,
              ToleranceToErrorRatio const& tolerance_to_error_ratio,
              Parameters const& parameters) const = 0;

  virtual void WriteToMessage(
      not_null<serialization::AdaptiveStepSizeIntegrator*> message) const = 0;
  static AdaptiveStepSizeIntegrator const& ReadFromMessage(
      serialization::AdaptiveStepSizeIntegrator const& message);

 protected:
  AdaptiveStepSizeIntegrator() = default;
};

template<typename Equation>
AdaptiveStepSizeIntegrator<Equation> const& ParseAdaptiveStepSizeIntegrator(
    std::string const& integrator_kind);

}  // namespace internal

using internal::AdaptiveStepSizeIntegrator;
using internal::FixedStepSizeIntegrator;
using internal::Integrator;
using internal::ParseAdaptiveStepSizeIntegrator;
using internal::ParseFixedStepSizeIntegrator;

}  // namespace _integrators
}  // namespace integrators
}  // namespace principia

namespace principia::integrators {
using namespace principia::integrators::_integrators;
}  // namespace principia::integrators

#include "integrators/integrators_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
