// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=ApproximationBenchmark  // NOLINT(whitespace/line_length)

#include <memory>
#include <random>
#include <vector>

#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/instant.hpp"
#include "numerics/approximation.hpp"
#include "numerics/polynomial_in_чебышёв_basis.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {

using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::numerics::_approximation;
using namespace principia::numerics::_polynomial_in_чебышёв_basis;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

constexpr int64_t number_of_angular_frequencies = 5;

template<int max_degree>
class ApproximationBenchmark : public benchmark::Fixture {
 protected:
  void SetUp(benchmark::State& state) override {
    std::mt19937_64 random(42);
    std::uniform_real_distribution<> coefficients(0.0, 10.0);
    std::uniform_real_distribution<> angular_frequency(0.0, 100.0);
    for (std::int64_t i = 0; i < number_of_angular_frequencies; ++i) {
      cos_coefficients_.push_back(coefficients(random) * Metre);
      sin_coefficients_.push_back(coefficients(random) * Metre);
      angular_frequencies_.push_back(angular_frequency(random) * Radian /
                                     Second);
    }
  }

  Length ComputeAltitude(Instant const& t) {
    Length altitude{};
    for (std::int64_t i = 0; i < number_of_angular_frequencies; ++i) {
      auto const& angular_frequency = angular_frequencies_[i];
      altitude += cos_coefficients_[i] * Cos(angular_frequency * (t - t0_)) +
                  sin_coefficients_[i] * Sin(angular_frequency * (t - t0_));
    }
    return altitude;
  }

  void ComputeAdaptiveЧебышёвPolynomialInterpolant(benchmark::State& state) {
    auto const altitude = [this](Instant const& t) {
      return ComputeAltitude(t);
    };
    SubdivisionPredicate<Length, Instant> const subdivide =
        [](auto const& interpolant, Length const& error_estimate) -> bool {
      return interpolant.MayHaveRealRoots(error_estimate);
    };

    std::vector<
        not_null<std::unique_ptr<PolynomialInЧебышёвBasis<Length, Instant>>>>
        interpolants;
    for (auto _ : state) {
      interpolants = AdaptiveЧебышёвPolynomialInterpolant<max_degree>(
          altitude,
          /*lower_bound=*/t0_,
          /*upper_bound=*/t0_ + 1000.0 * Second,
          /*max_error=*/1.0 * Metre,
          subdivide,
          /*error_estimate=*/nullptr);
      benchmark::DoNotOptimize(interpolants);
    }
  }

  Instant const t0_;
  std::vector<Length> cos_coefficients_;
  std::vector<Length> sin_coefficients_;
  std::vector<AngularFrequency> angular_frequencies_;
};

BENCHMARK_TEMPLATE_F(ApproximationBenchmark,
                     AdaptiveЧебышёвPolynomialInterpolant4,
                     4)(benchmark::State& state) {
  ComputeAdaptiveЧебышёвPolynomialInterpolant(state);
}
BENCHMARK_TEMPLATE_F(ApproximationBenchmark,
                     AdaptiveЧебышёвPolynomialInterpolant8,
                     8)(benchmark::State& state) {
  ComputeAdaptiveЧебышёвPolynomialInterpolant(state);
}
BENCHMARK_TEMPLATE_F(ApproximationBenchmark,
                     AdaptiveЧебышёвPolynomialInterpolant16,
                     16)(benchmark::State& state) {
  ComputeAdaptiveЧебышёвPolynomialInterpolant(state);
}
BENCHMARK_TEMPLATE_F(ApproximationBenchmark,
                     AdaptiveЧебышёвPolynomialInterpolant32,
                     32)(benchmark::State& state) {
  ComputeAdaptiveЧебышёвPolynomialInterpolant(state);
}

}  // namespace numerics
}  // namespace principia
