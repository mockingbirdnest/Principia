#pragma once

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/oblate_body.hpp
#include "quantities/named_quantities.hpp""
#include "quantities/tuples.hpp"

namespace principia {
namespace physics {
namespace internal_geopotential {

using base::not_null;
using geometry::Displacement;
using geometry::Instant;
using geometry::Vector;
using quantities::Acceleration;
using quantities::Angle;
using quantities::Derivatives;
using quantities::Exponentiation;
using quantities::GravitationalParameter;
using quantities::Infinity;
using quantities::Inverse;
using quantities::Length;
using quantities::Quotient;
using quantities::Square;

// Specification of the damping of a spherical harmonic, acting as a radial
// multiplier on the potential:
//   V_damped = σ(|r|) V(r).
struct HarmonicDamping final {
  explicit HarmonicDamping() = default;
  explicit HarmonicDamping(Length const& inner_threshold);
  
  // Above this threshold, the contribution to the potential from this
  // harmonic is 0, i.e., σ = 0.
  Length outer_threshold = Infinity<Length>();
  // Below this threshold, the contribution to the potential from this
  // harmonic is undamped, σ = 1.
  // inner_threshold = outer_threshold / 3.
  Length inner_threshold = Infinity<Length>();
  // For r in [outer_threshold, inner_threshold], σ is a polynomial with the
  // following coefficients in monomial basis.
  // The constant term is always 0, and is thus ignored in the evaluation.
  Derivatives<double, Length, 4> sigmoid_coefficients;

  // Sets σℜ_over_r and grad_σℜ according to σ as defined by |*this|.
  template<typename Frame>
  void ComputeDampedRadialQuantities(
      Length const& r_norm,
      Square<Length> const& r²,
      Vector<double, Frame> const& r_normalized,
      Inverse<Square<Length>> const& ℜ_over_r,
      Inverse<Square<Length>> const& ℜʹ,
      Inverse<Square<Length>>& σℜ_over_r,
      Vector<Inverse<Square<Length>>, Frame>& grad_σℜ) const;
};

// Representation of the geopotential model of an oblate body.
template<typename Frame>
class Geopotential {
 public:
  // Spherical harmonics will not be damped if they contribution to the radial
  // force exceeds |tolerance| times the central force.
  Geopotential(not_null<OblateBody<Frame> const*> body,
               double tolerance);

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  SphericalHarmonicsAcceleration(
      Instant const& t,
      Displacement<Frame> const& r,
      Square<Length> const& r²,
      Exponentiation<Length, -3> const& one_over_r³) const;

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  GeneralSphericalHarmonicsAcceleration(
      Instant const& t,
      Displacement<Frame> const& r,
      Length const& r_norm,
      Square<Length> const& r²,
      Exponentiation<Length, -3> const& one_over_r³) const;

  std::vector<HarmonicDamping> const& degree_damping() const;
  HarmonicDamping const& tesseral_damping() const;
  int first_tesseral_degree() const;

 private:
  // The frame of the surface of the celestial.
  struct SurfaceFrame;
  static const Vector<double, SurfaceFrame> x_;
  static const Vector<double, SurfaceFrame> y_;

  // This is the type that we return, so better have a name for it.
  using ReducedAcceleration = Quotient<Acceleration, GravitationalParameter>;

  // List of reduced accelerations computed for all degrees or orders.
  template<int size>
  using ReducedAccelerations =
      std::array<Vector<ReducedAcceleration, Frame>, size>;

  using UnitVector = Vector<double, Frame>;

  // Holds precomputed data for one evaluation of the acceleration.
  template<int size>
  struct Precomputations;

  // Helper templates for iterating over the degrees/orders of the geopotential.
  template<int size, int degree, int order>
  struct DegreeNOrderM;
  template<int size, int degree, typename>
  struct DegreeNAllOrders;
  template<typename>
  struct AllDegrees;

  // If z is a unit vector along the axis of rotation, and r a vector from the
  // center of |body_| to some point in space, the acceleration computed here
  // is:
  //
  //   -(J₂ / (μ ‖r‖⁵)) (3 z (r.z) + r (3 - 15 (r.z)² / ‖r‖²) / 2)
  //
  // Where ‖r‖ is the norm of r and r.z is the inner product.  It is the
  // additional acceleration exerted by the oblateness of |body| on a point at
  // position r.  J₂, J̃₂ and J̄₂ are normally positive and C̃₂₀ and C̄₂₀ negative
  // because the planets are oblate, not prolate.  Note that this follows IERS
  // Technical Note 36 and it differs from
  // https://en.wikipedia.org/wiki/Geopotential_model which seems to want J̃₂ to
  // be negative.
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Degree2ZonalAcceleration(UnitVector const& axis,
                           Displacement<Frame> const& r,
                           Exponentiation<Length, -2> const& one_over_r²,
                           Exponentiation<Length, -3> const& one_over_r³) const;

  not_null<OblateBody<Frame> const*> body_;

  // The contribution from the harmonics of degree n is damped by
  // degree_damping_[n].
  // degree_damping_[0] and degree_damping_[1] have infinite thresholds, and are
  // not used (this class does not compute the central force and disregards
  // degree 1, which is equivalent to a translation of the centre of mass).
  std::vector<HarmonicDamping> degree_damping_;

  // The contribution of low-degree tesseral (including sectoral) harmonics is
  // damped by |tesseral_damping_|.
  HarmonicDamping tesseral_damping_;

  // |first_tesseral_degree_| is the integer n such that
  //   degree_damping_[n-1].outer_threshold > tesseral_damping_.outer_threshold
  // and
  //   tesseral_damping_.outer_threshold >= degree_damping_[n].outer_threshold,
  // or is 0 if, for all n,
  //   degree_threshold_[n].outer_threshold <= tesseral_damping_.outer_threshold,
  // or |degree_damping_.size()| if, for all n,
  //   tesseral_damping_.outer_threshold < degree_threshold_[n].outer_threshold,
  // Tesseral (including sectoral) harmonics of degree less than
  // |first_tesseral_degree_| are damped by |tesseral_damping_|, instead
  // of their respective |degree_damping_|.
  int first_tesseral_degree_;
};

}  // namespace internal_geopotential

using internal_geopotential::Geopotential;

}  // namespace physics
}  // namespace principia

#include "physics/geopotential_body.hpp"
