#pragma once

#include<vector>

#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

using principia::quantities::Force;
using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::Momentum;
using principia::quantities::SIUnit;
using principia::quantities::Speed;
using principia::quantities::Stiffness;
using principia::quantities::Time;

namespace principia {
namespace testing_utilities {

inline void ComputeHarmonicOscillatorForce(Time const& t,
                                           std::vector<Length> const& q,
                                           std::vector<Force>* result) {
  (*result)[0] = -q[0] * SIUnit<Stiffness>();
}

inline void ComputeHarmonicOscillatorVelocity(std::vector<Momentum> const& p,
                                              std::vector<Speed>* result) {
  (*result)[0] = p[0] / SIUnit<Mass>();
}

}  // namespace testing_utilities
}  // namespace principia
