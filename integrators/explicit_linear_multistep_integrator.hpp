// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
//  parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_EXPLICIT_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_EXPLICIT_LINEAR_MULTISTEP_INTEGRATOR_HPP_

#include <list>
#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "integrators/starter.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace internal_explicit_linear_multistep_integrator {

using base::is_instance_of_v;
using base::not_null;
using numerics::DoublePrecision;
using numerics::FixedVector;

// This implementation follows [HW10], but restricts itself to explicit
// integrators.
template<typename Method, typename ODE_>
class ExplicitLinearMultistepIntegrator
    : public FixedStepSizeIntegrator<ODE_> {
 public:
  using ODE = ODE_;
  static_assert(
      is_instance_of_v<ExplicitFirstOrderOrdinaryDifferentialEquation, ODE>);
  using AppendState = typename Integrator<ODE>::AppendState;
  using IndependentVariable = typename ODE::IndependentVariable;
  using IndependentVariableDifference =
      typename ODE::IndependentVariableDifference;

  static constexpr int order = Method::order;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    absl::Status Solve(IndependentVariable const& s_final) override;
    ExplicitLinearMultistepIntegrator const& integrator() const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

   private:
    using DependentVariableDerivatives =
        typename ODE::DependentVariableDerivatives;
    using DoubleDependentVariables =
        decltype(std::declval<typename ODE::State>().y);

    // The data for a previous step of the integration.
    struct Step final {
      DoublePrecision<IndependentVariable> s;
      DoubleDependentVariables y;
      DependentVariableDerivatives yʹ;
    };

    class Starter : public integrators::Starter<ODE, Step, /*steps=*/order> {
     protected:
      using integrators::Starter<ODE, Step, order>::Starter;

      void FillStepFromState(ODE const& equation,
                             typename ODE::State const& state,
                             Step& step) const override;
    };

    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state,
             IndependentVariableDifference const& step,
             ExplicitLinearMultistepIntegrator const& integrator);

    Starter starter_;
    ExplicitLinearMultistepIntegrator const& integrator_;
    friend class ExplicitLinearMultistepIntegrator;
  };

  explicit ExplicitLinearMultistepIntegrator(
      FixedStepSizeIntegrator<ODE> const& startup_integrator);

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      InitialValueProblem<ODE> const& problem,
      AppendState const& append_state,
      IndependentVariableDifference const& step) const override;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const override;

 private:
  // See [HW10], equation (7).
  static constexpr auto α_ = Method::α;
  static constexpr auto β_numerator_ = Method::β_numerator;
  static constexpr auto β_denominator_ = Method::β_denominator;

  FixedStepSizeIntegrator<ODE> const& startup_integrator_;
};

}  // namespace internal_explicit_linear_multistep_integrator

template<typename Method, typename Position>
internal_explicit_linear_multistep_integrator::
    ExplicitLinearMultistepIntegrator<Method, Position> const&
ExplicitLinearMultistepIntegrator();

}  // namespace integrators
}  // namespace principia

#include "explicit_linear_multistep_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_EXPLICIT_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
