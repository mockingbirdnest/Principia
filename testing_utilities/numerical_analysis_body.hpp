#pragma once

#include <vector>

#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using quantities::Force;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;
using quantities::Momentum;
using quantities::Pow;
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

inline void ComputeKeplerAcceleration(
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>* const result) {
  auto const r_squared = q[0] * q[0] + q[1] * q[1];
  auto const minus_μ_over_r_cubed =
      -SIUnit<GravitationalParameter>() * Sqrt(r_squared) /
          (r_squared * r_squared);
  (*result)[0] = q[0] * minus_μ_over_r_cubed;
  (*result)[1] = q[1] * minus_μ_over_r_cubed;
}

}  // namespace testing_utilities
}  // namespace principia
