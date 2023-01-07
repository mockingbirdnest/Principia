// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_INTEGRATOR_HPP_

#include <functional>
#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "numerics/fixed_arrays.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_runge_kutta_integrator {

using base::is_instance_of_v;
using base::not_null;
using numerics::FixedStrictlyLowerTriangularMatrix;
using numerics::FixedVector;
using quantities::Variation;

// This class solves ordinary differential equations of the form q″ = f(q, t)
// using an embedded Runge-Kutta method.  We follow the standard conventions for
// the coefficients, i.e.,
//   c for the nodes;
//   a for the Runge-Kutta matrix;
//   b̂ for the weights of the high-order method;
//   b for the weights of the low-order method;
// See [DP86] for an example.

// In the implementation, we follow [DP807] in calling the results of the
// right-hand-side evaluations fᵢ (this quantity is not named in [DP86]).
// The order of the template parameters follow the notation of [DP86], whose
// RKq(p)s[F]X has higher order q, lower order p, comprises s stages, and has
// the first-same-as-last property.

template<typename Method, typename ODE_>
class EmbeddedExplicitRungeKuttaIntegrator
    : public AdaptiveStepSizeIntegrator<ODE_> {
 public:
  using ODE = ODE_;
  static_assert(
      is_instance_of_v<ExplicitFirstOrderOrdinaryDifferentialEquation, ODE>);
  using typename Integrator<ODE>::AppendState;
  using typename AdaptiveStepSizeIntegrator<ODE>::Parameters;
  using typename AdaptiveStepSizeIntegrator<ODE>::ToleranceToErrorRatio;

  static constexpr auto higher_order = Method::higher_order;
  static constexpr auto lower_order = Method::lower_order;
  static constexpr auto first_same_as_last = Method::first_same_as_last;

  EmbeddedExplicitRungeKuttaIntegrator();

  EmbeddedExplicitRungeKuttaIntegrator(
      EmbeddedExplicitRungeKuttaIntegrator const&) = delete;
  EmbeddedExplicitRungeKuttaIntegrator(
      EmbeddedExplicitRungeKuttaIntegrator&&) = delete;
  EmbeddedExplicitRungeKuttaIntegrator& operator=(
      EmbeddedExplicitRungeKuttaIntegrator const&) = delete;
  EmbeddedExplicitRungeKuttaIntegrator& operator=(
      EmbeddedExplicitRungeKuttaIntegrator&&) = delete;

  class Instance : public AdaptiveStepSizeIntegrator<ODE>::Instance {
   public:
    absl::Status Solve(
        typename ODE::IndependentVariable const& s_final) override;
    EmbeddedExplicitRungeKuttaIntegrator const& integrator()
        const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
#if 0
    template<typename... DV = DependentVariable...,
             typename = std::enable_if_t<base::is_serializable_v<DV...>>>
    static not_null<std::unique_ptr<Instance>> ReadFromMessage(
        serialization::
            EmbeddedExplicitRungeKuttaNystromIntegratorInstance const&
                extension,
        InitialValueProblem<ODE> const& problem,
        AppendState const& append_state,
        ToleranceToErrorRatio const& tolerance_to_error_ratio,
        Parameters const& parameters,
        typename ODE::IndependentVariableDifference const& step,
        bool first_use,
        EmbeddedExplicitRungeKuttaIntegrator const& integrator);
#endif

   private:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state,
             ToleranceToErrorRatio const& tolerance_to_error_ratio,
             Parameters const& parameters,
             typename ODE::IndependentVariableDifference const& step,
             bool first_use,
             EmbeddedExplicitRungeKuttaIntegrator const& integrator);

    EmbeddedExplicitRungeKuttaIntegrator const& integrator_;
    friend class EmbeddedExplicitRungeKuttaIntegrator;
  };

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      InitialValueProblem<ODE> const& problem,
      AppendState const& append_state,
      ToleranceToErrorRatio const& tolerance_to_error_ratio,
      Parameters const& parameters) const override;

  void WriteToMessage(
      not_null<serialization::AdaptiveStepSizeIntegrator*> message)
      const override;

 private:
  static constexpr auto stages_ = Method::stages;
  static constexpr auto c_ = Method::c;
  static constexpr auto a_ = Method::a;
  static constexpr auto b̂_ = Method::b̂;
  static constexpr auto b_ = Method::b;
};

}  // namespace internal_embedded_explicit_runge_kutta_integrator

template<typename Method, typename ODE>
internal_embedded_explicit_runge_kutta_integrator::
    EmbeddedExplicitRungeKuttaIntegrator<Method, ODE> const&
EmbeddedExplicitRungeKuttaIntegrator();

}  // namespace integrators
}  // namespace principia

#include "integrators/embedded_explicit_runge_kutta_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
