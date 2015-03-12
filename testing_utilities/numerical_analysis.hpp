#pragma once

#include<vector>

#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using quantities::Acceleration;
using quantities::Force;
using quantities::Length;
using quantities::Momentum;
using quantities::Speed;
using quantities::Time;

namespace testing_utilities {

// Right-hand sides for various differential equations frequently used to test
// the properties of integrators.  Logically the result should be a not_null<>
// but this costs 20% in the benchmark of the harmonic oscillator.

// The one-dimensional unit harmonic oscillator,
// q' = p,  |ComputeHarmonicOscillatorVelocity|,
// p' = -q, |ComputeHarmonicOscillatorForce|.

void ComputeHarmonicOscillatorForce(Time const& t,
                                    std::vector<Length> const& q,
                                    std::vector<Force>* const result);

void ComputeHarmonicOscillatorVelocity(
    std::vector<Momentum> const& p,
    std::vector<Speed>* const result);

// The Runge-Kutta-Nyström formulation.
void ComputeHarmonicOscillatorAcceleration(
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>* const result);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/numerical_analysis_body.hpp"
