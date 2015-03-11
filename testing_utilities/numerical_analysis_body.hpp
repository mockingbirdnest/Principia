#pragma once

#include<vector>

#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using quantities::Force;
using quantities::Length;
using quantities::Mass;
using quantities::Momentum;
using quantities::SIUnit;
using quantities::Speed;
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

}  // namespace testing_utilities
}  // namespace principia
