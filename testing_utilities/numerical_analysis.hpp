#pragma once

#include<vector>

// Right-hand sides for various differential equations frequently used to test
// the properties of integrators.

namespace principia {
namespace testing_utilities {

// The one-dimensional unit harmonic oscillator,
// q' = p,  |ComputeHarmonicOscillatorVelocity|,
// p' = -q, |ComputeHarmonicOscillatorForce|.

void ComputeHarmonicOscillatorVelocity(std::vector<double> const& p,
                                       std::vector<double>* result);

void ComputeHarmonicOscillatorForce(double const t,
                                    std::vector<double> const& q,
                                    std::vector<double>* result);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/numerical_analysis_body.hpp"
