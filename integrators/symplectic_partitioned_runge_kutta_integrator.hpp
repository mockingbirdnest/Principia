
#pragma once

#include <type_traits>

#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_runge_kutta_nyström_integrator {

using numerics::FixedVector;

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
class SymplecticPartitionedRungeKuttaIntegrator {
 public:
  SymplecticPartitionedRungeKuttaIntegrator();
};

}  // namespace internal_symplectic_runge_kutta_nyström_integrator

using internal_symplectic_runge_kutta_nyström_integrator::
    SymplecticPartitionedRungeKuttaIntegrator;

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
