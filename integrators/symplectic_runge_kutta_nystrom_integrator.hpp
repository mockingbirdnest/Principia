
#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "integrators/symplectic_integrator.hpp"

namespace principia {

using base::not_null;

namespace integrators {

// The integrator defined in this class integrates Hamilton's equations for
// Hamiltonians of the form
//   H(q, p, t) = T(q, p) + V(q, t),
// with T quadratic in p, i.e.,
//  T(q, p) = ½ M(q)⁻¹(p, p)
// where M(q)⁻¹ is bilinear.  The flows of T and V must be known.
// A special treatment is provided for hamiltonians where M⁻¹ is independent of
// q and ∂V/∂qᵢ(q, t) is known.  In those cases the differential equation can be
// reformulated as q" = M⁻¹ F(q, t), where Fᵢ = ∂V/∂qᵢ.

// NOTE(egg): only the latter is implemented at this time.

// TODO(egg): a blurb on non-autonomy, it's slightly trickier than in the SPRK
// case.

// We follow the convention of McLachlan & Atela, calling the position nodes
// aᵢ and the momentum nodes bᵢ.  The following notations appear in the
// litterature:
//   (a, b) in McLachlan & Atela, as well as Candy & Rozmus;
//   (d, c) in Ruth, Yoshida, as well as Forest & Ruth;
//   (B, b) in Sofroniou & Spaletta.
// Moreover, we follow the convention of Sofroniou & Spaletta in calling the
// weights used for the time argument of the force computation cᵢ, with
// c₁ = 0, cᵢ = cᵢ₋₁ + aᵢ₋₁ for i > 1.

enum VanishingCoefficients {
  kNone,
  kFirstBVanishes,
  kLastAVanishes,
};

class SRKNIntegrator : public SymplecticIntegrator {
 protected:
  enum VanishingCoefficients {
   kNone,
   kFirstBVanishes,
   kLastAVanishes,
 };

 public:

  // The functors |evolve_kinetic| and |evolve_potential| respectively compute
  // exp(h{·, T})(q, p) given (q, p), and exp(h{·, V})(q, p) given (q, p, t).
  template<typename Position, typename Momentum,
           typename KineticFlow
           typename PotentialFlow>
  void SolveQuadraticKineticEnergy(
      KineticFlow evolve_kinetic,
      PotentialFlow evolve_potential,
      Parameters<Position, Momentum> const& parameters,
      not_null<Solution<Position, Momentum>*> const solution) const;

  // The functor |compute_acceleration| computes M⁻¹ F(q, t).
  template<typename Position, typename Velocity,
           typename RightHandSideComputation>
  void SolveTrivialKineticEnergyIncrement(
      RightHandSideComputation compute_acceleration,
      Parameters<Position, Velocity> const& parameters,
      not_null<Solution<Position, Velocity>*> const solution) const;

 protected:
  struct FirstSameAsLast {
    double first;
    double last;
  };

  VanishingCoefficients vanishing_coefficients_;
  // Null if, and only if, |vanishing_coefficients_ == kNone|.
  // If |vanishing_coefficients_ == kFirstBVanishes|, this contains the first
  // and last a coefficients.
  // If |vanishing_coefficients_ == FirstAVanishes|, this contains the first
  // and last b coefficients.
  // These are used to desynchronize and synchronize first-same-as-last
  // integrators.
  std::unique_ptr<FirstSameAsLast> first_same_as_last_

  int const stages_;

  // The position and momentum nodes.
  // If |vanishing_coefficients_ == kFirstBVanishes|, we do not store the first
  // b coefficient, and the last entry of |a_| is the sum of the first and last
  // a coefficients.
  std::vector<double> a_;
  std::vector<double> b_;

  // The weights.  Note that the first c coefficient is not stored if
  // |vanishing_coefficients_ == kFirstBVanishes|.
  std::vector<double> c_;
 private:
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
