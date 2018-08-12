#pragma once

#include "physics/geopotential.hpp"

namespace principia {
namespace physics {
namespace internal_geopotential {

template<typename Frame>
Geopotential<Frame>::Geopotential(not_null<OblateBody<Frame>* const> body)
    : body_(body) {}

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::SphericalHarmonicsAcceleration(
    Instant const& t,
    Displacement<Frame> const& r) {
  return 
}

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::Order2ZonalAcceleration(Displacement<Frame> const& r) {}

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
