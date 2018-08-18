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

 private:
  // The frame of the surface of the celestial.
  struct SurfaceFrame;
  static const Vector<double, SurfaceFrame> x_;
  static const Vector<double, SurfaceFrame> y_;

  using UnitVector = Vector<double, Frame>;

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

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Degree2SectoralAcceleration(
      UnitVector const& reference,
      UnitVector const& bireference,
      Displacement<Frame> const& r,
      Exponentiation<Length, -2> const& one_over_r²,
      Exponentiation<Length, -3> const& one_over_r³) const;

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Degree3ZonalAcceleration(UnitVector const& axis,
                           Displacement<Frame> const& r,
                           Square<Length> const& r²,
                           Exponentiation<Length, -2> const& one_over_r²,
                           Exponentiation<Length, -3> const& one_over_r³) const;

  not_null<OblateBody<Frame> const*> const body_;
};

}  // namespace internal_geopotential

using internal_geopotential::Geopotential;

}  // namespace physics
}  // namespace principia

#include "physics/geopotential_body.hpp"
