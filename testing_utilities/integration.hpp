#pragma once

#include <vector>

#include "absl/status/status.h"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "physics/massive_body.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

#include <tuple>

namespace principia {
namespace testing_utilities {
namespace _integration {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::physics::_massive_body;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// Right-hand sides for various differential equations frequently used to test
// the properties of integrators.

// The Runge-Kutta-Nyström formulation
//   qʺ = -q k / m.
absl::Status ComputeHarmonicOscillatorAcceleration1D(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>& result,
    int* evaluations);

template<typename Frame>
absl::Status ComputeHarmonicOscillatorAcceleration3D(
    Instant const& t,
    std::vector<Position<Frame>> const& q,
    std::vector<Vector<Acceleration, Frame>>& result,
    int* evaluations);

// The Runge-Kutta formulation
//   qʹ = v
//   vʹ = -q k / m.
absl::Status ComputeHarmonicOscillatorDerivatives1D(
    Instant const& t,
    std::tuple<Length, Speed> const& state,
    std::tuple<Speed, Acceleration>& result,
    int* evaluations);

// The Kepler problem with unit gravitational parameter, where the
// two-dimensional configuration space is the separation between the bodies, in
// the Runge-Kutta-Nyström formulation
//   qʺ = -q μ / |q|³,
// where μ = 1 m³ s⁻².
absl::Status ComputeKeplerAcceleration(Instant const& t,
                                       std::vector<Length> const& q,
                                       std::vector<Acceleration>& result,
                                       int* evaluations);

// The right-hand side of the Чебышёв differential equation, with the
// independent variable scaled so that the interval [-1, 1] maps to
// [J2000 - 1 s, J2000 + 1 s].
// ч, чʹ, and чʺ must have size 1.
template<int degree>
absl::Status ComputeЧебышёвPolynomialSecondDerivative(
    Instant const& t,
    std::vector<double> const& ч,
    std::vector<Variation<double>> const& чʹ,
    std::vector<Variation<Variation<double>>>& чʺ,
    int* evaluations);

// The right-hand side of the Legendre differential equation, with the
// independent variable scaled so that the interval [-1, 1] maps to
// [J2000 - 1 s, J2000 + 1 s].
// p, pʹ, and pʺ must have size 1.
template<int degree>
absl::Status ComputeLegendrePolynomialSecondDerivative(
    Instant const& t,
    std::vector<double> const& p,
    std::vector<Variation<double>> const& pʹ,
    std::vector<Variation<Variation<double>>>& pʺ,
    int* evaluations);

}  // namespace internal

using internal::ComputeЧебышёвPolynomialSecondDerivative;
using internal::ComputeHarmonicOscillatorAcceleration1D;
using internal::ComputeHarmonicOscillatorAcceleration3D;
using internal::ComputeHarmonicOscillatorDerivatives1D;
using internal::ComputeKeplerAcceleration;
using internal::ComputeLegendrePolynomialSecondDerivative;

}  // namespace _integration
}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/integration_body.hpp"
