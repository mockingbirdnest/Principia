
// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)
#define PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)

#include <functional>
#include <vector>

#include "base/not_null.hpp"
#include "base/status.hpp"
#include "numerics/fixed_arrays.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_runge_kutta_nyström_integrator {

using base::not_null;
using base::Status;
using geometry::Instant;
using numerics::FixedStrictlyLowerTriangularMatrix;
using numerics::FixedVector;
using quantities::Time;
using quantities::Variation;

// This class solves ordinary differential equations of the form q″ = f(q, t)
// using an embedded Runge-Kutta-Nyström method.  We follow the standard
// conventions for the coefficients, i.e.,
//   c for the nodes;
//   a for the Runge-Kutta matrix;
//   b̂ for the position weights of the high-order method;
//   b̂′ for the velocity weights of the high-order method;
//   b for the position weights of the low-order method;
//   b′ for the velocity weights of the low-order method.
// See Dormand, El-Mikkawy and Prince (1986),
// Families of Runge-Kutta-Nyström formulae, for an example.
// In the implementation, we follow Dormand, El-Mikkawy and Prince in calling
// the results of the right-hand-side evaluations gᵢ.  The notations kᵢ or fᵢ
// also appear in the litterature.
// Since we are interested in physical applications, we call the solution q and
// its derivative v, rather than the more common y and y′ found in the
// literature on Runge-Kutta-Nyström methods.
// The order of the template parameters follow the notation of Dormand and
// Prince, whose RKNq(p)sF has higher order q, lower order p, comprises
// s stages, and has the first-same-as-last property.

template<typename Method, typename Position>
class EmbeddedExplicitRungeKuttaNyströmIntegrator
    : public AdaptiveStepSizeIntegrator<
                 SpecialSecondOrderDifferentialEquation<Position>> {
 public:
  using ODE = SpecialSecondOrderDifferentialEquation<Position>;
  using typename Integrator<ODE>::AppendState;
  using typename AdaptiveStepSizeIntegrator<ODE>::Parameters;
  using typename AdaptiveStepSizeIntegrator<ODE>::ToleranceToErrorRatio;

  static constexpr auto higher_order = Method::higher_order;
  static constexpr auto lower_order = Method::lower_order;
  static constexpr auto first_same_as_last = Method::first_same_as_last;

  EmbeddedExplicitRungeKuttaNyströmIntegrator();

  EmbeddedExplicitRungeKuttaNyströmIntegrator(
      EmbeddedExplicitRungeKuttaNyströmIntegrator const&) = delete;
  EmbeddedExplicitRungeKuttaNyströmIntegrator(
      EmbeddedExplicitRungeKuttaNyströmIntegrator&&) = delete;
  EmbeddedExplicitRungeKuttaNyströmIntegrator& operator=(
      EmbeddedExplicitRungeKuttaNyströmIntegrator const&) = delete;
  EmbeddedExplicitRungeKuttaNyströmIntegrator& operator=(
      EmbeddedExplicitRungeKuttaNyströmIntegrator&&) = delete;

  class Instance : public AdaptiveStepSizeIntegrator<ODE>::Instance {
   public:
    Status Solve(Instant const& t_final) override;
    EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator()
        const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;

   private:
    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state,
             ToleranceToErrorRatio const& tolerance_to_error_ratio,
             Parameters const& adaptive_step_size,
             EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator);

    EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator_;
    friend class EmbeddedExplicitRungeKuttaNyströmIntegrator;
  };

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      ToleranceToErrorRatio const& tolerance_to_error_ratio,
      Parameters const& parameters) const override;

 private:
  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>
  ReadFromMessage(
      serialization::AdaptiveStepSizeIntegratorInstance const& message,
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      ToleranceToErrorRatio const& tolerance_to_error_ratio,
      Parameters const& parameters) const override;

  static constexpr auto stages_ = Method::stages;
  static constexpr auto c_ = Method::c;
  static constexpr auto a_ = Method::a;
  static constexpr auto b_hat_ = Method::b_hat;
  static constexpr auto b_prime_hat_ = Method::b_prime_hat;
  static constexpr auto b_ = Method::b;
  static constexpr auto b_prime_ = Method::b_prime;
};

}  // namespace internal_embedded_explicit_runge_kutta_nyström_integrator

template<typename Method, typename Position>
internal_embedded_explicit_runge_kutta_nyström_integrator::
    EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position> const&
EmbeddedExplicitRungeKuttaNyströmIntegrator();

}  // namespace integrators
}  // namespace principia

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
