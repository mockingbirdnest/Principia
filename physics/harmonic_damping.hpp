#pragma once

#include "geometry/grassmann.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace _harmonic_damping {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_quantities;

// Specification of the damping of a spherical harmonic, acting as a radial
// multiplier on the potential:
//   V_damped = σ(‖r‖) V(r).
class HarmonicDamping final {
 public:
  HarmonicDamping() = default;
  explicit HarmonicDamping(Length const& inner_threshold);

  // Above this threshold, the contribution to the potential from this
  // harmonic is 0, i.e., σ = 0.
  Length const& outer_threshold() const;
  // Below this threshold, the contribution to the potential from this
  // harmonic is undamped, σ = 1.
  // This class depends on the invariant: outer_threshold = 3 * inner_threshold.
  Length const& inner_threshold() const;

  // Sets σℜ_over_r and grad_σℜ according to σ as defined by `*this`.
  template<typename Frame>
  void ComputeDampedRadialQuantities(
      Length const& r_norm,
      Square<Length> const& r²,
      Vector<double, Frame> const& r_normalized,
      Inverse<Square<Length>> const& ℜ_over_r,
      Inverse<Square<Length>> const& ℜʹ,
      Inverse<Square<Length>>& σℜ_over_r,
      Vector<Inverse<Square<Length>>, Frame>& grad_σℜ) const;

  // Same as above, but only computes the quantities needed for the potential.
  void ComputeDampedRadialQuantities(Length const& r_norm,
                                     Square<Length> const& r²,
                                     Inverse<Square<Length>> const& ℜ_over_r,
                                     Inverse<Square<Length>>& σℜ_over_r) const;

 private:
  Length outer_threshold_ = Infinity<Length>;
  Length inner_threshold_ = Infinity<Length>;

  // For r in [outer_threshold, inner_threshold], σ is a polynomial with the
  // following coefficients in monomial basis.
  // The constant term is always 0, and is thus ignored in the evaluation.
  // TODO(phl): We don't use an evaluator; we use a custom evaluation that
  // ignores the constant term instead.  See #1922.
  PolynomialInMonomialBasis<double, Length, 3>::Coefficients
      sigmoid_coefficients_;
};
}  // namespace internal

using internal::HarmonicDamping;

}  // namespace _harmonic_damping
}  // namespace physics
}  // namespace principia

#include "physics/harmonic_damping_body.hpp"
