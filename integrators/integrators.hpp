#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#define PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_

#include <functional>
#include <limits>

#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace internal_integrators {

using base::Error;
using base::not_null;
using base::Status;
using geometry::Instant;
using numerics::DoublePrecision;
using quantities::Time;

// A base class for integrators.
template<typename ODE_>
class Integrator {
 public:
  using ODE = ODE_;
  using AppendState =
      std::function<void(typename ODE::SystemState const& state)>;

  // An object for holding the integrator state during the integration of a
  // problem.
  class Instance {
   public:
    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state);
    virtual ~Instance() = default;

    // The subclass must document the time passed to the last call to
    // |append_state|.
    virtual Status Solve(Instant const& t_final) = 0;

    // The last instant integrated by this instance.
    DoublePrecision<Instant> const& time() const;

    // The last state integrated by this instance.
    typename ODE::SystemState const& state() const;

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
    typename ODE::SystemState current_state_;
    AppendState const append_state_;
  };

  virtual ~Integrator() = default;
};

// An integrator using a fixed step size.
template<typename ODE_>
class FixedStepSizeIntegrator : public Integrator<ODE_> {
 public:
  using ODE = ODE_;

  // The last call to |append_state| has a |state.time.value| equal to the
  // unique |Instant| of the form |problem.t_final + n * step| in
  // ]t_final - step, t_final].  |append_state| will be called with
  // |state.time.values|s at intervals differing from |step| by at most one ULP.
  class Instance : public Integrator<ODE>::Instance {
   public:
    // The integrator corresponding to this instance.
    virtual FixedStepSizeIntegrator const& integrator() const = 0;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
    static not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
    ReadFromMessage(serialization::IntegratorInstance const& message,
                    ODE equation,
                    AppendState const& append_state);

   protected:
    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state,
             Time const& step);

    Time const step_;
  };

  // The factory function for |Instance|, above.  It ensures that the instance
  // has a back-pointer to its integrator.
  virtual not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
  NewInstance(IntegrationProblem<ODE> const& problem,
              typename Integrator<ODE>::AppendState const& append_state,
              Time const& step) const = 0;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const;
  static FixedStepSizeIntegrator const& ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message);

 protected:
  // For convenience, deserialization is an instance member of the |Integrator|,
  // not a static member of the |Instance|.  Which makes sense if you see
  // |Integrator| as a factory for |Instance|.
  virtual not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
  ReadFromMessage(
      serialization::FixedStepSizeIntegratorInstance const& message,
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const = 0;

  explicit FixedStepSizeIntegrator(
      serialization::FixedStepSizeIntegrator::Kind kind);

 private:
  serialization::FixedStepSizeIntegrator::Kind const kind_;
};

// An integrator using an adaptive step size.
template<typename ODE_>
class AdaptiveStepSizeIntegrator : public Integrator<ODE_> {
 public:
  using ODE = ODE_;

  struct Parameters final {
    using ToleranceToErrorRatio =
        std::function<
            double(Time const& current_step_size,
                   typename ODE::SystemStateError const& error)>;
    // The first time step tried by the integrator. It must have the same sign
    // as |problem.t_final - initial_state.time.value|.
    Time first_time_step;
    // This number must be in ]0, 1[.  Higher values increase the chance of step
    // rejection, lower values yield smaller steps.
    double safety_factor;
    // This functor is called at each step, with the |current_step_size| used by
    // the integrator and the estimated |error| on that step.  It returns the
    // ratio of a tolerance to some norm of the error.  The step is recomputed
    // with a smaller step size if the result is less than 1, and accepted
    // otherwise.
    // In both cases, the new step size is chosen so as to try and make the
    // result of the next call to |tolerance_to_error_ratio| close to
    // |safety_factor|.
    ToleranceToErrorRatio tolerance_to_error_ratio;
    // Integration will stop after |*max_steps| even if it has not reached
    // |t_final|.
    std::int64_t max_steps = std::numeric_limits<std::int64_t>::max();

    void WriteToMessage(
        not_null<serialization::AdaptiveStepSizeIntegratorInstance::
                     Parameters*> const message) const;
    static Parameters ReadFromMessage(
        serialization::AdaptiveStepSizeIntegratorInstance::Parameters const&
            message);
  };

  // The last call to |append_state| will have |state.time.value == t_final|.
  class Instance : public Integrator<ODE>::Instance {
   public:
    // The integrator corresponding to this instance.
    virtual AdaptiveStepSizeIntegrator const& integrator() const = 0;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
    static not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
    ReadFromMessage(serialization::IntegratorInstance const& message,
                    ODE equation,
                    AppendState const& append_state,
                    typename Parameters::ToleranceToErrorRatio const&
                        tolerance_to_error_ratio);

   protected:
    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state,
             Parameters const& parameters);

    Parameters const parameters_;
  };

  // The factory function for |Instance|, above.  It ensures that the instance
  // has a back-pointer to its integrator.
  virtual not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
  NewInstance(IntegrationProblem<ODE> const& problem,
              typename Integrator<ODE>::AppendState const& append_state,
              Parameters const& parameters) const = 0;

  void WriteToMessage(
      not_null<serialization::AdaptiveStepSizeIntegrator*> message) const;
  static AdaptiveStepSizeIntegrator const& ReadFromMessage(
      serialization::AdaptiveStepSizeIntegrator const& message);

 protected:
  // For convenience, deserialization is an instance member of the |Integrator|,
  // not a static member of the |Instance|.  Which makes sense if you see
  // |Integrator| as a factory for |Instance|.
  virtual not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
  ReadFromMessage(
      serialization::AdaptiveStepSizeIntegratorInstance const& message,
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      Parameters const& parameters) const = 0;

  explicit AdaptiveStepSizeIntegrator(
      serialization::AdaptiveStepSizeIntegrator::Kind kind);

 private:
  serialization::AdaptiveStepSizeIntegrator::Kind const kind_;
};

}  // namespace internal_integrators

using internal_integrators::AdaptiveStepSizeIntegrator;
using internal_integrators::FixedStepSizeIntegrator;
using internal_integrators::Integrator;

}  // namespace integrators
}  // namespace principia

#include "integrators/integrators_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
