
#pragma once

#include <vector>

#include "testing_utilities/integration.hpp"

#include "astronomy/epoch.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_integration {

using astronomy::J2000;
using geometry::Displacement;
using geometry::InnerProduct;
using quantities::Exponentiation;
using quantities::Force;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;
using quantities::Momentum;
using quantities::Pow;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Square;
using quantities::Stiffness;
using quantities::Time;
using quantities::si::Second;

inline Status ComputeHarmonicOscillatorAcceleration1D(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>& result,
    int* const evaluations) {
  result[0] = -q[0] * (SIUnit<Stiffness>() / SIUnit<Mass>());
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return Status::OK;
}

template<typename Frame>
Status ComputeHarmonicOscillatorAcceleration3D(
    Instant const& t,
    std::vector<Position<Frame>> const& q,
    std::vector<Vector<Acceleration, Frame>>& result,
    int* const evaluations) {
  result[0] = (Frame::origin - q[0]) * (SIUnit<Stiffness>() / SIUnit<Mass>());
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return Status::OK;
}

inline Status ComputeKeplerAcceleration(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>& result,
    int* evaluations) {
  auto const r² = q[0] * q[0] + q[1] * q[1];
  auto const minus_μ_over_r³ =
      -SIUnit<GravitationalParameter>() * Sqrt(r²) / (r² * r²);
  result[0] = q[0] * minus_μ_over_r³;
  result[1] = q[1] * minus_μ_over_r³;
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return Status::OK;
}


template<int degree>
Status ComputeЧебышёвPolynomialSecondDerivative(
    Instant const& t,
    std::vector<double> const& ч,
    std::vector<Variation<double>> const& чʹ,
    std::vector<Variation<Variation<double>>>& чʺ,
    int* evaluations) {
  constexpr int n² = degree * degree;
  constexpr auto s² = Pow<2>(Second);
  Time const x = (t - J2000);
  auto const x² = x * x;
  чʺ[0] = (x * чʹ[0] - n² * ч[0]) / (1 * s² - x²);
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return Status::OK;
}

template<int degree>
Status ComputeLegendrePolynomialSecondDerivative(
    Instant const& t,
    std::vector<double> const& p,
    std::vector<Variation<double>> const& pʹ,
    std::vector<Variation<Variation<double>>>& pʺ,
    int* evaluations) {
  constexpr int n = degree;
  constexpr auto s² = Pow<2>(Second);
  Time const x = (t - J2000);
  auto const x² = x * x;
  pʺ[0] = (2 * x * pʹ[0] - n * (n + 1) * p[0]) / (1 * s² - x²);
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return Status::OK;
}

}  // namespace internal_integration
}  // namespace testing_utilities
}  // namespace principia
