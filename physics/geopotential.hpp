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
using quantities::GravitationalParameter;
using quantities::Quotient;

template<typename Frame>
class Geopotential {
 public:
  explicit Geopotential(not_null<OblateBody<Frame>* const> body);

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  SphericalHarmonicsAcceleration(Instant const& t,
                                 Displacement<Frame> const& r);

 private:
  using UnitVector = Vector<double, Frame>;

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Order2ZonalAcceleration(UnitVector const& i,
                          UnitVector const& j,
                          UnitVector const& k,
                          Displacement<Frame> const& r);

  not_null<OblateBody<Frame>* const> const body_;
};

}  // namespace internal_geopotential

using internal_geopotential::Geopotential;

}  // namespace physics
}  // namespace principia

#include "physics/geopotential_body.hpp"
