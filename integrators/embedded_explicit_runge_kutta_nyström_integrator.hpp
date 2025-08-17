// The files containing the tree of of child classes of `Integrator` must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)
#define PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)

#include <functional>
#include <memory>

#include "absl/status/status.h"
#include "base/concepts.hpp"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/instant.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace _embedded_explicit_runge_kutta_nyström_integrator {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_instant;
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_quantities;

// This class solves ordinary differential equations of the form q″ = f(q, t)
// using an embedded Runge-Kutta-Nyström method.  We follow the standard
// conventions for the coefficients, i.e.,
//   c for the nodes;
//   a for the position Runge-Kutta matrix;
//   b̂ for the position weights of the high-order method;
//   b̂′ for the velocity weights of the high-order method;
//   b for the position weights of the low-order method;
//   b′ for the velocity weights of the low-order method.
// See [DEP87a], for an example.
// Note that other notations exist for the weights: [ER03] and [Som93] call the
// velocity weights d instead of b′, and [ACM06] call the position and velocity
// weights β and b instead of b and b′.  Further, [Dor96] uses ā for the
// Runge-Kutta matrix, and b̄ and b for the position and velocity weights.

// In the implementation, we follow Dormand, El-Mikkawy and Prince in calling
// the results of the right-hand-side evaluations gᵢ.  The notations kᵢ or fᵢ
// also appear in the literature.
// Since we are interested in physical applications, we call the solution q and
// its derivative v, rather than the more common y and y′ found in the
// literature on Runge-Kutta-Nyström methods.
// The order of the template parameters follow the notation of [DEP87a], which
// is a generalization of that of [DP86]: the method RKNq(p)s[F]X has higher
// order q, lower order p, comprises s stages, and has the first-same-as-last
// property.

template<typename Method, typename ODE_>
class EmbeddedExplicitRungeKuttaNyströmIntegrator
    : public AdaptiveStepSizeIntegrator<ODE_> {
 public:
  using ODE = ODE_;
  static_assert(is_instance_of_v<SpecialSecondOrderDifferentialEquation, ODE>);
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
    absl::Status Solve(Instant const& t_final) override;
    EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator()
        const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
    static not_null<std::unique_ptr<Instance>> ReadFromMessage(
        serialization::
            EmbeddedExplicitRungeKuttaNystromIntegratorInstance const&
                extension,
        InitialValueProblem<ODE> const& problem,
        AppendState const& append_state,
        ToleranceToErrorRatio const& tolerance_to_error_ratio,
        Parameters const& parameters,
        Time const& time_step,
        bool first_use,
        EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator)
      requires serializable<typename ODE_::DependentVariable>;

   private:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state,
             ToleranceToErrorRatio const& tolerance_to_error_ratio,
             Parameters const& parameters,
             Time const& time_step,
             bool first_use,
             EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator);

    EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator_;
    friend class EmbeddedExplicitRungeKuttaNyströmIntegrator;
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
  static constexpr auto b̂ʹ_ = Method::b̂ʹ;
  static constexpr auto b_ = Method::b;
  static constexpr auto bʹ_ = Method::bʹ;
};

}  // namespace internal

template<typename Method, typename Position>
internal::EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position> const&
EmbeddedExplicitRungeKuttaNyströmIntegrator();

}  // namespace _embedded_explicit_runge_kutta_nyström_integrator
}  // namespace integrators
}  // namespace principia

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
