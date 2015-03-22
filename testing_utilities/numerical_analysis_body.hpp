#pragma once

#include<vector>

#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using geometry::Displacement;
using geometry::InnerProduct;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using quantities::Area;
using quantities::Force;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;
using quantities::Momentum;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Stiffness;
using quantities::Time;

namespace testing_utilities {

inline void ComputeHarmonicOscillatorForce(
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Force>* const result) {
  (*result)[0] = -q[0] * SIUnit<Stiffness>();
}

inline void ComputeHarmonicOscillatorVelocity(
    std::vector<Momentum> const& p,
    std::vector<Speed>* const result) {
  (*result)[0] = p[0] / SIUnit<Mass>();
}

inline void ComputeHarmonicOscillatorAcceleration(
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>* const result) {
  (*result)[0] = -q[0] * (SIUnit<Stiffness>() / SIUnit<Mass>());
}

template<typename Frame>
inline void ComputeKeplerAcceleration(
    Time const& t,
    std::vector<Position<Frame>> const& q,
    std::vector<Vector<Acceleration, Frame>>* const result) {
  Displacement<Frame> const r = q[1] - q[0];
  Area const r_squared = InnerProduct(r, r);
  auto const μ_over_r_cubed =
      (0.5 * SIUnit<GravitationalParameter>() * Sqrt(r_squared)) /
          (r_squared * r_squared);
  (*result)[0] = r * μ_over_r_cubed;
  (*result)[1] = -(*result)[0];
}

}  // namespace testing_utilities
}  // namespace principia
