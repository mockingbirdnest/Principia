#pragma once

#include<vector>

#include "base/not_null.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using base::not_null;
using quantities::Force;
using quantities::Length;
using quantities::Momentum;
using quantities::Speed;
using quantities::Time;

namespace testing_utilities {

// Right-hand sides for various differential equations frequently used to test
// the properties of integrators.

// The one-dimensional unit harmonic oscillator,
// q' = p,  |ComputeHarmonicOscillatorVelocity|,
// p' = -q, |ComputeHarmonicOscillatorForce|.

void ComputeHarmonicOscillatorForce(Time const& t,
                                    std::vector<Length> const& q,
                                    not_null<std::vector<Force>*> const result);

void ComputeHarmonicOscillatorVelocity(
    std::vector<Momentum> const& p,
    not_null<std::vector<Speed>*> const result);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/numerical_analysis_body.hpp"
