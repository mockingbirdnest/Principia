
// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
#include "integrators/ordinary_differential_equations.hpp"
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

template<typename Position, int higher_order, int lower_order, int stages,
         bool first_same_as_last>
class EmbeddedExplicitRungeKuttaNyströmIntegrator
    : public AdaptiveStepSizeIntegrator<
                 SpecialSecondOrderDifferentialEquation<Position>> {
 public:
  using ODE = SpecialSecondOrderDifferentialEquation<Position>;

  EmbeddedExplicitRungeKuttaNyströmIntegrator(
      serialization::AdaptiveStepSizeIntegrator::Kind kind,
      FixedVector<double, stages> const& c,
      FixedStrictlyLowerTriangularMatrix<double, stages> const& a,
      FixedVector<double, stages> const& b_hat,
      FixedVector<double, stages> const& b_prime_hat,
      FixedVector<double, stages> const& b,
      FixedVector<double, stages> const& b_prime);

  EmbeddedExplicitRungeKuttaNyströmIntegrator() = delete;
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
    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;

   private:
    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state,
             AdaptiveStepSize<ODE> const& adaptive_step_size,
             EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator);

    EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator_;
    friend class EmbeddedExplicitRungeKuttaNyströmIntegrator;
  };

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
    IntegrationProblem<ODE> const& problem,
    typename Integrator<ODE>::AppendState const& append_state,
    AdaptiveStepSize<ODE> const& adaptive_step_size) const override;

  // The actual integration, applied to |instance|.
  Status Solve(Instant const& t_final, Instance& instance) const;

 protected:
  FixedVector<double, stages> const c_;
  FixedStrictlyLowerTriangularMatrix<double, stages> const a_;
  FixedVector<double, stages> const b_hat_;
  FixedVector<double, stages> const b_prime_hat_;
  FixedVector<double, stages> const b_;
  FixedVector<double, stages> const b_prime_;
};

// Coefficients from Dormand, El-Mikkawy and Prince (1986),
// Families of Runge-Kutta-Nyström formulae, table 3 (the RK4(3)4FM).
// Minimizes the 4th order truncation error.
template<typename Position>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Position,
                                            /*higher_order=*/4,
                                            /*lower_order=*/3,
                                            /*stages=*/4,
                                            /*first_same_as_last=*/true> const&
DormandElMikkawyPrince1986RKN434FM();

}  // namespace internal_embedded_explicit_runge_kutta_nyström_integrator

using internal_embedded_explicit_runge_kutta_nyström_integrator::
    DormandElMikkawyPrince1986RKN434FM;
using internal_embedded_explicit_runge_kutta_nyström_integrator::
    EmbeddedExplicitRungeKuttaNyströmIntegrator;

}  // namespace integrators
}  // namespace principia

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_EMBEDDED_EXPLICIT_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_  // NOLINT(whitespace/line_length)
#endif  // PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
