
#pragma once

#include "numerics/fixed_arrays.hpp"
#include "symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_partitioned_runge_kutta_integrator {

using numerics::FixedVector;

// A symplectic partitioned Runge-Kutta integrator.  Does not subclass
// |Integrator| yet; used to generate (less general)
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
// |SymplecticRungeKuttaNyströmIntegrator|; since A and B are interchangeable
// for a |SymplecticPartitionedRungeKuttaIntegrator|, this may be done by either
// making B the "force operator" and A the "velocity operator", corresponding to
//   [B, [B, [B, A]]] = 0,
// or by making A the "force operator" and B the "velocity operator",
// corresponding to
//   [A, [A, [A, B]]] = 0.
// If the method is |first_same_as_last|, the former yields a |BAB|
// |CompositionMethod|, and the latter yields an |ABA| |CompositionMethod|.
// If the method is not |first_same_as_last|, both yield a |BA|
// |CompositionMethod|.
// NOTE(egg): The |SymplecticRungeKuttaNyströmIntegrator| thus constructed will
// serialize as a |DUMMY| and probably break in all sorts of hilarious ways if
// deserialized.
// TODO(egg): Make them serializable/deserializable.  We need to prevent
// combinatorial explosion.
template<typename Position, typename Momentum,
         int order_, int evaluations_, bool time_reversible_,
         bool first_same_as_last_>
class SymplecticPartitionedRungeKuttaIntegrator {
  static constexpr int stages_ = first_same_as_last_ ? evaluations_ + 1
                                                     : evaluations_;
 public:
  SymplecticPartitionedRungeKuttaIntegrator(
      FixedVector<double, stages_> const& a,
      FixedVector<double, stages_> const& b);

  static constexpr int order = order_;
  static constexpr int evaluations = evaluations_;
  static constexpr bool time_reversible = time_reversible_;
  static constexpr bool first_same_as_last = first_same_as_last_;

  SymplecticRungeKuttaNyströmIntegrator<
      Position,
      order,
      time_reversible,
      evaluations,
      first_same_as_last ? BAB : BA> const& BForceMethod() const;

  SymplecticRungeKuttaNyströmIntegrator<
      Position,
      order,
      time_reversible,
      evaluations,
      first_same_as_last ? BAB : BA> const& AForceMethod() const;

 private:
  FixedVector<double, stages_> const a_;
  FixedVector<double, stages_> const b_;

  // The Runge-Kutta-Nyström methods are stored here, so that we can use them
  // by const-reference as we do for the others.  Since |*this| should be a
  // static object, This similarly obviates questions of lifetime.
  std::unique_ptr<SymplecticRungeKuttaNyströmIntegrator<
      Position,
      order,
      time_reversible,
      evaluations,
      first_same_as_last_ ? BAB : BA>> b_force_method_;
  std::unique_ptr<SymplecticRungeKuttaNyströmIntegrator<Position,
      order,
      time_reversible,
      evaluations,
      first_same_as_last_ ? ABA : BA>> a_force_method_;
};

}  // namespace internal_symplectic_partitioned_runge_kutta_integrator

using internal_symplectic_partitioned_runge_kutta_integrator::
    SymplecticPartitionedRungeKuttaIntegrator;

}  // namespace integrators
}  // namespace principia
