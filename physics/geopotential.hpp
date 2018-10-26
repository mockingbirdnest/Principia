#pragma once

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/oblate_body.hpp"
#include "quantities/named_quantities.hpp"

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
using quantities::Length;
using quantities::Quotient;
using quantities::Square;

// Representation of the geopotential model of an oblate body.
template<typename Frame>
class Geopotential {
 public:
  explicit Geopotential(not_null<OblateBody<Frame> const*> body);

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

  // Beyond |degree_threshold_[n]|, the contribution from the harmonics with
  // degree ≥n is ignored.  Between |degree_threshold / 3| and
  // |degree_threshold|, a sigmoid is applied to the relevant terms of the
  // potential. In particular, beyond |degree_threshold_[2]|,
  // |GeneralSphericalHarmonicsAcceleration| returns 0.
  // |degree_threshold_[0]| and |degree_threshold_[1]| are unused, and are
  // infinite.
  // The values in |degree_threshold_| are in non-strictly decreasing order.
  std::vector<Length> degree_threshold_;

  // The coefficients of the sigmoid, in monomial basis.  A custom evaluation is
  // performed to take advantage of precomputations.  The constant term must be
  // 0: it is not used in the evaluation.
  std::vector<Derivatives<double, Length, 4>> degree_sigmoid_coefficients_;

  // Beyond this threshold, the contribution from the tesseral harmonics
  // (including the sectoral harmonics) is 0.  Between
  // |tesseral_threshold_ / 3| and |tesseral_threshold_|, a sigmoid is applied.
  Length tesseral_threshold_;

  // The coefficients of the sigmoid, in monomial basis.  A custom evaluation is
  // performed to take advantage of precomputations.  The constant term must be
  // 0: it is not used in the evaluation.
  Derivatives<double, Length, 4> tesseral_sigmoid_coefficients_;

  // |first_tesseral_degree_| is the integer n such that
  // |degree_threshold_[n-1] >= tesseral_threshold_ > degree_threshold_[n]|,
  // or is |degree_threshold_.size()| if
  // |tesseral_threshold_ <= degree_threshold_[n]| for all n.
  // Tesseral (including sectoral) harmonics of degree less than
  // |first_tesseral_degree_| come into effect at |tesseral_threshold_|, instead
  // of their respective |degree_threshold_|.
  int first_tesseral_degree_;
};

}  // namespace internal_geopotential

using internal_geopotential::Geopotential;

}  // namespace physics
}  // namespace principia

#include "physics/geopotential_body.hpp"
