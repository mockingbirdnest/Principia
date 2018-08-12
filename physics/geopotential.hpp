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

template<typename Frame>
class Geopotential {
 public:
  explicit Geopotential(not_null<OblateBody<Frame> const*> body);

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  SphericalHarmonicsAcceleration(Instant const& t,
                                 Displacement<Frame> const& r,
                                 Square<Length> const& r²,
                                 Exponentiation<Length, -3> const& one_over_r³);

 private:
  using UnitVector = Vector<double, Frame>;

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Order2ZonalAcceleration(UnitVector const& axis,
                          Displacement<Frame> const& r,
                          Exponentiation<Length, -2> const& one_over_r²,
                          Exponentiation<Length, -3> const& one_over_r³);

  not_null<OblateBody<Frame> const*> const body_;
};

}  // namespace internal_geopotential

using internal_geopotential::Geopotential;

}  // namespace physics
}  // namespace principia

#include "physics/geopotential_body.hpp"
