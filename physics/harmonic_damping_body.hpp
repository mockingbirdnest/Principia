#pragma once

#include "physics/harmonic_damping.hpp"

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_harmonic_damping {

using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

// The notation in this file follows documentation/Geopotential.pdf.

inline HarmonicDamping::HarmonicDamping(Length const& inner_threshold)
    : outer_threshold_(inner_threshold * 3),
      inner_threshold_(inner_threshold),
      sigmoid_coefficients_{0,
                            9 / (4 * inner_threshold),
                            -3 / (2 * Pow<2>(inner_threshold)),
                            1 / (4 * Pow<3>(inner_threshold))} {}

inline Length const& HarmonicDamping::outer_threshold() const {
  return outer_threshold_;
}

inline Length const& HarmonicDamping::inner_threshold() const {
  return inner_threshold_;
}

template<typename Frame>
void HarmonicDamping::ComputeDampedRadialQuantities(
      Length const& r_norm,
      Square<Length> const& r²,
      Vector<double, Frame> const& r_normalized,
      Inverse<Square<Length>> const& ℜ_over_r,
      Inverse<Square<Length>> const& ℜʹ,
      Inverse<Square<Length>>& σℜ_over_r,
      Vector<Inverse<Square<Length>>, Frame>& grad_σℜ) const {
  Length const& s0 = inner_threshold_;
  if (r_norm <= s0) {
    // Below the inner threshold, σ = 1.
    σℜ_over_r = ℜ_over_r;
    grad_σℜ = ℜʹ * r_normalized;
  } else {
    auto const& c = sigmoid_coefficients_;
    Derivative<double, Length> const c1 = std::get<1>(c);
    Derivative<double, Length, 2> const c2 = std::get<2>(c);
    Derivative<double, Length, 3> const c3 = std::get<3>(c);
    auto const r³ = r² * r_norm;
    double const c3r³ = c3 * r³;
    double const c2r² = c2 * r²;
    double const c1r = c1 * r_norm;
    double const σ = c3r³ + c2r² + c1r;
    double const σʹr = 3 * c3r³ + 2 * c2r² + c1r;

    σℜ_over_r = σ * ℜ_over_r;
    // Writing this as σ′ℜ + ℜ′σ rather than ℜ∇σ + σ∇ℜ turns some vector
    // operations into scalar ones.
    grad_σℜ = (σʹr * ℜ_over_r + ℜʹ * σ) * r_normalized;
  }
}

inline void HarmonicDamping::ComputeDampedRadialQuantities(
    Length const& r_norm,
    Square<Length> const& r²,
    Inverse<Square<Length>> const& ℜ_over_r,
    Inverse<Square<Length>>& σℜ_over_r) const {
  Length const& s0 = inner_threshold_;
  if (r_norm <= s0) {
    // Below the inner threshold, σ = 1.
    σℜ_over_r = ℜ_over_r;
  } else {
    auto const& c = sigmoid_coefficients_;
    Derivative<double, Length> const c1 = std::get<1>(c);
    Derivative<double, Length, 2> const c2 = std::get<2>(c);
    Derivative<double, Length, 3> const c3 = std::get<3>(c);
    auto const r³ = r² * r_norm;
    double const c3r³ = c3 * r³;
    double const c2r² = c2 * r²;
    double const c1r = c1 * r_norm;
    double const σ = c3r³ + c2r² + c1r;

    σℜ_over_r = σ * ℜ_over_r;
  }
}

}  // namespace internal_harmonic_damping
}  // namespace physics
}  // namespace principia
