#pragma once

#include "base/not_null.hpp"
#include "integrators/integrators.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/double_precision.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_starter {

using base::not_null;
using numerics::DoublePrecision;
using quantities::Difference;

// A helper object for starting a linear multistep integrator.  |ODE| is the
// equation being integrated.  |Step| is an object holding the data produced by
// executing a step of the integrator and used for subsequent steps.  The
// function |FillStepFromState| must fill this object from an |ODE::State|.
// |order| is the order of the integrator being started.
template<typename ODE, typename Step, int order>
class Starter {
 public:
  // |startup_step_divisor| specifies how many steps of the |startup_integrator|
  // should be produced for a step of the main integrator.  |instance| is the
  // instance of the main integrator.
  Starter(FixedStepSizeIntegrator<ODE> const& startup_integrator,
          std::int64_t startup_step_divisor,
          not_null<typename FixedStepSizeIntegrator<ODE>::Instance*> const
              instance);

  // Performs the startup integration, i.e., computes enough states to either
  // reach |s_final| or to reach a point where |instance.previous_steps_| has
  // |order - 1| elements.  During startup |instance.current_state_| is
  // updated more frequently than once every |instance.step_|.
  void StartupSolve(typename ODE::IndependentVariable const& s_final);

 protected:
  // Must fill |step| from |state| and the right-hand side of the ODE.  |Step|
  // must contain all the information needed to compute subsequent steps of the
  // integrator.
  virtual void FillStepFromState(
      typename ODE::RightHandSideComputation const& rhs,
      typename ODE::State const& state,
      Step& step) = 0;

 private:
  FixedStepSizeIntegrator<ODE> const& startup_integrator_;
  std::int64_t const startup_step_divisor_;
  not_null<typename FixedStepSizeIntegrator<ODE>::Instance*> const
      instance_;

  int startup_step_index_ = 0;
  std::list<Step> previous_steps_;  // At most |order_| elements.
};

}  // namespace internal_starter

using internal_starter::Starter;

}  // namespace integrators
}  // namespace principia

#include "integrators/starter_body.hpp"
