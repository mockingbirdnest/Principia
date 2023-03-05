// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_GENERALIZED_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)
#define PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_GENERALIZED_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)

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
namespace _embedded_explicit_generalized_runge_kutta_nyström_integrator {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_named_quantities;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// This class solves ordinary differential equations of the form
//   q″ = f(q, q′, t)
// using an embedded generalized Runge-Kutta-Nyström (RKNG) method.  We follow
// the following conventions for the coefficients:
//   c for the nodes;
//   a for the position Runge-Kutta matrix;
//   a′ for the velocity Runge-Kutta matrix;
//   b̂ for the position weights of the high-order method;
//   b̂′ for the velocity weights of the high-order method;
//   b for the position weights of the low-order method;
//   b′ for the velocity weights of the low-order method.
// This is similar to the notation used, e.g., by [Fin87].  This notation is
// consistent with the one we use for Runge-Kutta-Nyström methods (which have no
// a′).  However, note that Fine uses ^ for the low order method, whereas we use
// it for the high order method, consistently with the notation used by Dormand,
// المكاوى, and Prince for RKN methods.

// Alternative notations use a and b for the velocity Runge Kutta matrix and
// weights: [Mur98] uses (c, α, a, β̂, b̂, β, b), and [Dor96] uses
// (c, ā, a, b̄̂, b̂, b̄, b).
// These alternative notations are more consistent with the usual notation for
// Runge-Kutta methods, as a Runge-Kutta method with nodes cᵢ, Runge-Kutta
// matrix aⁱⱼ, and weights bᵢ, applied to the 1st order ODE
// (q, q′)′ = (q′, f(q, q′, t)), is equivalent to a Runge-Kutta-Nyström
// Generalized method with, using Einstein summation convention,
//   nodes cᵢ;
//   position Runge-Kutta matrix aⁱₖaᵏⱼ;
//   velocity Runge-Kutta matrix aⁱⱼ;
//   position weights bₖaᵏᵢ;
//   velocity weights bᵢ.

template<typename Method, typename ODE_>
class EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator
    : public AdaptiveStepSizeIntegrator<ODE_> {
 public:
  using ODE = ODE_;
  static_assert(
      is_instance_of_v<ExplicitSecondOrderOrdinaryDifferentialEquation, ODE>);
  using typename Integrator<ODE>::AppendState;
  using typename AdaptiveStepSizeIntegrator<ODE>::Parameters;
  using typename AdaptiveStepSizeIntegrator<ODE>::ToleranceToErrorRatio;

  static constexpr auto higher_order = Method::higher_order;
  static constexpr auto lower_order = Method::lower_order;
  static constexpr auto first_same_as_last = Method::first_same_as_last;

  EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator();

  EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator(
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator const&) = delete;
  EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator(
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator&&) = delete;
  EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator& operator=(
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator const&) = delete;
  EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator& operator=(
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator&&) = delete;

  class Instance : public AdaptiveStepSizeIntegrator<ODE>::Instance {
   public:
    absl::Status Solve(Instant const& t_final) override;
    EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator const& integrator()
        const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
    template<typename DV = typename ODE::DependentVariable,
             typename = std::enable_if_t<base::is_serializable_v<DV>>>
    static not_null<std::unique_ptr<Instance>> ReadFromMessage(
        serialization::
        EmbeddedExplicitGeneralizedRungeKuttaNystromIntegratorInstance const&
                extension,
        InitialValueProblem<ODE> const& problem,
        AppendState const& append_state,
        ToleranceToErrorRatio const& tolerance_to_error_ratio,
        Parameters const& parameters,
        Time const& time_step,
        bool first_use,
        EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator const&
            integrator);

   private:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state,
             ToleranceToErrorRatio const& tolerance_to_error_ratio,
             Parameters const& parameters,
             Time const& time_step,
             bool first_use,
             EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator const&
                 integrator);

    EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator const& integrator_;
    friend class EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator;
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
  static constexpr auto aʹ_ = Method::aʹ;
  static constexpr auto b̂_ = Method::b̂;
  static constexpr auto b̂ʹ_ = Method::b̂ʹ;
  static constexpr auto b_ = Method::b;
  static constexpr auto bʹ_ = Method::bʹ;
};

}  // namespace internal

template<typename Method, typename ODE>
internal_embedded_explicit_generalized_runge_kutta_nyström_integrator::
    EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, ODE> const&
EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator();

}  // namespace _embedded_explicit_generalized_runge_kutta_nyström_integrator
}  // namespace integrators
}  // namespace principia

namespace principia::integrators {
using namespace principia::integrators::_embedded_explicit_generalized_runge_kutta_nyström_integrator;
}  // namespace principia::integrators

#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator_body.hpp"  // NOLINT(whitespace/line_length)

#endif  // PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_GENERALIZED_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
