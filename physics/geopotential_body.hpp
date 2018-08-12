#pragma once

#include "physics/geopotential.hpp"

#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_geopotential {

using geometry::InnerProduct;
using quantities::Length;

template<typename Frame>
Geopotential<Frame>::Geopotential(not_null<OblateBody<Frame> const*> body)
    : body_(body) {}

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::SphericalHarmonicsAcceleration(
    Instant const& t,
    Displacement<Frame> const& r,
    Square<Length> const& r²,
    Exponentiation<Length, -3> const& one_over_r³) {
  Exponentiation<Length, -2> const one_over_r² = 1 / r²;
  Vector<double, Frame> const& axis = body_->polar_axis();
  return Order2ZonalAcceleration(axis, r, one_over_r², one_over_r³);
}

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::Order2ZonalAcceleration(
    UnitVector const& axis,
    Displacement<Frame> const& r,
    Exponentiation<Length, -2> const& one_over_r²,
    Exponentiation<Length, -3> const& one_over_r³) {
  Length const r_axis_projection = InnerProduct(axis, r);
  auto const j2_over_r_fifth = body_->j2_over_μ() * one_over_r³ * one_over_r²;
  Vector<Quotient<Acceleration,
                  GravitationalParameter>, Frame> const axis_effect =
      (-3 * j2_over_r_fifth * r_axis_projection) * axis;
  Vector<Quotient<Acceleration,
                  GravitationalParameter>, Frame> const radial_effect =
      (j2_over_r_fifth *
           (-1.5 +
            7.5 * r_axis_projection * r_axis_projection * one_over_r²)) * r;
  return axis_effect + radial_effect;
}

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
