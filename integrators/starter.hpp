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

// |order| is the order of the integrator being started.
template<typename ODE, int order>
class Starter {
 public:
  using typename ODE::IndependentVariable;
  using typename ODE::DependentVariable;

  Starter(FixedStepSizeIntegrator<ODE> const& startup_integrator,
          not_null<typename FixedStepSizeIntegrator<ODE>::Instance*> const
              instance);

  // Performs the startup integration, i.e., computes enough states to either
  // reach |s_final| or to reach a point where |instance.previous_steps_| has
  // |order - 1| elements.  During startup |instance.current_state_| is
  // updated more frequently than once every |instance.step_|.
  void StartupSolve(IndependentVariable const& s_final);

 private:
  struct Step final {
    std::vector<DoublePrecision<Difference<DependentVariable>>> displacements;
    std::vector<typename ODE::Acceleration> accelerations;
    DoublePrecision<IndependentVariable> s;

    void WriteToMessage(
        not_null<
            serialization::SymmetricLinearMultistepIntegratorInstance::Step*>
            message) const;
    template<typename P = Position,
             typename = std::enable_if_t<base::is_serializable_v<P>>>
    static Step ReadFromMessage(
        serialization::SymmetricLinearMultistepIntegratorInstance::Step const&
            message);
  };

  static void FillStepFromState(ODE const& equation,
                                typename ODE::State const& state,
                                Step& step);

  FixedStepSizeIntegrator<ODE> const& startup_integrator_;
  not_null<typename FixedStepSizeIntegrator<ODE>::Instance*> const instance_;

  int startup_step_index_ = 0;
  std::list<Step> previous_steps_;  // At most |order_| elements.
};

}  // namespace internal_starter

using internal_starter::Starter;

}  // namespace integrators
}  // namespace principia

#include "integrators/starter_body.hpp"
