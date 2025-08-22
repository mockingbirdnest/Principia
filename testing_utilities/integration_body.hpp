#pragma once

#include "testing_utilities/integration.hpp"

#include <tuple>
#include <vector>

#include "astronomy/epoch.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace testing_utilities {
namespace _integration {
namespace internal {

using namespace principia::astronomy::_epoch;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;

inline absl::Status ComputeHarmonicOscillatorAcceleration1D(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>& result,
    int* const evaluations) {
  result[0] = -q[0] * (si::Unit<Stiffness> / si::Unit<Mass>);
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return absl::OkStatus();
}

template<typename Frame>
absl::Status ComputeHarmonicOscillatorAcceleration3D(
    Instant const& t,
    std::vector<Position<Frame>> const& q,
    std::vector<Vector<Acceleration, Frame>>& result,
    int* const evaluations) {
  result[0] = (Frame::origin - q[0]) * (si::Unit<Stiffness> / si::Unit<Mass>);
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return absl::OkStatus();
}

inline absl::Status ComputeHarmonicOscillatorDerivatives1D(
    Instant const& t,
    std::tuple<Length, Speed> const& state,
    std::tuple<Speed, Acceleration>& result,
    int* const evaluations) {
  auto const& [q, v] = state;
  auto& [qʹ, vʹ] = result;
  qʹ = v;
  vʹ = -q * (si::Unit<Stiffness> / si::Unit<Mass>);
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return absl::OkStatus();
}

template<typename Frame>
absl::Status ComputeKeplerAcceleration(
    Instant const& t,
    std::vector<Position<Frame>> const& q,
    std::vector<Vector<Acceleration, Frame>>& result,
    int* evaluations) {
  auto const r = q[0] - Frame::origin;
  auto const r³ = Pow<3>(r.Norm());
  result[0] = -r * si::Unit<GravitationalParameter> / r³;
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return absl::OkStatus();
}


template<int degree>
absl::Status ComputeЧебышёвPolynomialSecondDerivative(
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
  return absl::OkStatus();
}

template<int degree>
absl::Status ComputeLegendrePolynomialSecondDerivative(
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
  return absl::OkStatus();
}

}  // namespace internal
}  // namespace _integration
}  // namespace testing_utilities
}  // namespace principia
