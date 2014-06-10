#pragma once

#include<vector>

#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

using principia::quantities::Force;
using principia::quantities::Length;
using principia::quantities::Momentum;
using principia::quantities::Speed;
using principia::quantities::Time;

// Right-hand sides for various differential equations frequently used to test
// the properties of integrators.

namespace principia {
namespace testing_utilities {

// The one-dimensional unit harmonic oscillator,
// q' = p,  |ComputeHarmonicOscillatorVelocity|,
// p' = -q, |ComputeHarmonicOscillatorForce|.

void ComputeHarmonicOscillatorForce(Time const& t,
                                    std::vector<Length> const& q,
                                    std::vector<Force>* result);

void ComputeHarmonicOscillatorVelocity(std::vector<Momentum> const& p,
                                       std::vector<Speed>* result);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/numerical_analysis_body.hpp"
