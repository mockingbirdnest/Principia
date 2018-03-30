
#pragma once

#include <type_traits>

#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/integrators.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_partitioned_runge_kutta_integrator {

using base::not_null;
using base::Status;
using geometry::Instant;
using quantities::Time;

// A symplectic partitioned Runge-Kutta integrator.  Does not subclass
// |Integrator|; used to generate (less general)
// |SymplecticRungeKuttaNyströmIntegrator|s.
// Represents a single-step method for the solution of
//   (q, p)′ = X(q, p, t), with X = A(q, p, t) + B(q, p, t).
// |Position| is the type of |q|, and |Momentum| is that of |p|.
// The step is the composition of evolutions
//   exp(aᵣ₋₁ h A) exp(bᵣ₋₁ h B) ... exp(a₀ h A) exp(b₀ h B);
// A and B are interchangeable.  If aᵣ₋₁ vanishes, this becomes
//   exp(bᵣ₋₁ h B) exp(aᵣ₋₂ h A) ... exp(a₀ h A) exp(b₀ h B).
// In that case, |first_same_as_last| is true.
// The equation solved by this integrator is a general case of that solved by
// a |SymplecticRungeKuttaNyströmIntegrator|, see equation (1) in the
// appropriate file.  This may therefore be turned into a
// |SymplecticRungeKuttaNyströmIntegrator|.
// In the |first_same_as_last| case, since A and B are interchangeable
// for a |SymplecticPartitionedRungeKuttaIntegrator|, the step
//   exp(bᵣ₋₁ h A) exp(aᵣ₋₂ h B) ... exp(a₀ h B) exp(b₀ h A)
// also provides an integrator.  Since for a
// |SymplecticRungeKuttaNyströmIntegrator| A and B are *not* interchangeable,
// because of the requirement that [B, [B, [B, A]]] = 0 in (1), there are two
// different |SymplecticRungeKuttaNyströmIntegrator| corresponding to a
// |first_same_as_last| |SymplecticPartitionedRungeKuttaIntegrator|: one is
// |ABA|, the other is |BAB|.
// NOTE(egg): The |SymplecticRungeKuttaNyströmIntegrator| thus constructed will
// serialize as a |DUMMY| and probably break in all sorts of hilarious ways if
// deserialized.
// TODO(egg): Make them serializable/deserializable.  We need to prevent
// combinatorial explosion.
template<typename Method, typename Position>
class SymplecticPartitionedRungeKuttaIntegrator
    : public FixedStepSizeIntegrator<
          DecomposableFirstOrderDifferentialEquation<Position>> {
 public:
  using ODE = DecomposableFirstOrderDifferentialEquation<Position>;
  using AppendState = typename Integrator<ODE>::AppendState;

  static constexpr auto time_reversible = Method::time_reversible;
  static constexpr auto first_same_as_last = Method::first_same_as_last;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    Status Solve(Instant const& t_final) override;
    SymplecticPartitionedRungeKuttaIntegrator const& integrator()
        const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;

   private:
    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state,
             Time const& step,
             SymplecticPartitionedRungeKuttaIntegrator const& integrator);

    SymplecticPartitionedRungeKuttaIntegrator const& integrator_;
    friend class SymplecticPartitionedRungeKuttaIntegrator;
  };

  SymplecticPartitionedRungeKuttaIntegrator();

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const override;

 private:
  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> ReadFromMessage(
      serialization::FixedStepSizeIntegratorInstance const& message,
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const override;

  static constexpr auto stages_ = Method::stages;
  static constexpr auto a_ = Method::a;
  static constexpr auto b_ = Method::b;
};

}  // namespace internal_symplectic_partitioned_runge_kutta_integrator

template<typename Method, typename Position>
internal_symplectic_partitioned_runge_kutta_integrator::
    SymplecticPartitionedRungeKuttaIntegrator<Method, Position> const&
SymplecticPartitionedRungeKuttaIntegrator();

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
