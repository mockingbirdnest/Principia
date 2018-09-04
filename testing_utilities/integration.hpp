
#pragma once

#include <vector>

#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/massive_body.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_integration {

using base::not_null;
using base::Status;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using physics::MassiveBody;
using quantities::Acceleration;
using quantities::Force;
using quantities::Length;
using quantities::Momentum;
using quantities::Speed;
using quantities::Time;
using quantities::Variation;

// Right-hand sides for various differential equations frequently used to test
// the properties of integrators.

// The Runge-Kutta-Nyström formulation
//   qʺ = -q k / m.
Status ComputeHarmonicOscillatorAcceleration1D(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>& result,
    int* evaluations);

template<typename Frame>
Status ComputeHarmonicOscillatorAcceleration3D(
    Instant const& t,
    std::vector<Position<Frame>> const& q,
    std::vector<Vector<Acceleration, Frame>>& result,
    int* evaluations);

// The Kepler problem with unit gravitational parameter, where the
// two-dimensional configuration space is the separation between the bodies, in
// the Runge-Kutta-Nyström formulation
//   q" = -q μ / |q|³,
// where μ = 1 m³ s⁻².
Status ComputeKeplerAcceleration(Instant const& t,
                                 std::vector<Length> const& q,
                                 std::vector<Acceleration>& result,
                                 int* evaluations);

// The right-hand side of the Чебышёв differential equation, with the
// independent variable scaled so that the interval [-1, 1] maps to
// [J2000 - 1 s, J2000 + 1 s].
// ч, чʹ, and чʺ must have size 1.
template<int degree>
Status ComputeЧебышёвPolynomialSecondDerivative(
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
Status ComputeLegendrePolynomialSecondDerivative(
    Instant const& t,
    std::vector<double> const& p,
    std::vector<Variation<double>> const& pʹ,
    std::vector<Variation<Variation<double>>>& pʺ,
    int* evaluations);

}  // namespace internal_integration

using internal_integration::ComputeЧебышёвPolynomialSecondDerivative;
using internal_integration::ComputeHarmonicOscillatorAcceleration1D;
using internal_integration::ComputeHarmonicOscillatorAcceleration3D;
using internal_integration::ComputeKeplerAcceleration;
using internal_integration::ComputeLegendrePolynomialSecondDerivative;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/integration_body.hpp"
