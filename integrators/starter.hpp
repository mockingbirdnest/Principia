#pragma once

#include "numerics/double_precision.hpp"

namespace principia {
namespace integrators {
namespace internal_starter {

using numerics::DoublePrecision;

template<typename IndependentVariable>
class Starter {
 public:
  // Performs the startup integration, i.e., computes enough states to either
  // reach |s_final| or to reach a point where |instance.previous_steps_| has
  // |order - 1| elements.  During startup |instance.current_state_| is
  // updated more frequently than once every |instance.step_|.
  void StartupSolve(IndependentVariable const& s_final);

 private:
  struct Step final {
    std::vector<DoublePrecision<typename ODE::Displacement>> displacements;
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

  static void FillStepFromSystemState(ODE const& equation,
                                      typename ODE::SystemState const& state,
                                      Step& step);
  int startup_step_index_ = 0;
  std::list<Step> previous_steps_;  // At most |order_| elements.
};

}  // namespace internal_starter

using internal_starter::Starter;

}  // namespace integrators
}  // namespace principia
