// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=EvaluatePolynomialInЧебышёвBasis  // NOLINT(whitespace/line_length)

#include <random>
#include <string>
#include <vector>

#include "astronomy/frames.hpp"
#include "benchmark/benchmark.h"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/space.hpp"
#include "numerics/polynomial_in_чебышёв_basis.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {

using namespace principia::astronomy::_frames;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_space;
using namespace principia::numerics::_polynomial_in_чебышёв_basis;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

namespace {
constexpr int evaluations_per_iteration = 1000;
}  // namespace

template<int degree>
void BM_EvaluatePolynomialInЧебышёвBasisDouble(benchmark::State& state) {
  std::mt19937_64 random(42);
  std::array<double, degree + 1> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients[i] = static_cast<double>(random());
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  PolynomialInЧебышёвBasis<double, Instant, degree> const series(
      coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  double result = 0.0;

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series(t);
      t += Δt;
    }
  }

  // This weird call to `SetLabel` has no effect except that it uses `result`
  // and therefore prevents the loop from being optimized away.
  state.SetLabel(std::to_string(result).substr(0, 0));
}

template<int degree>
void BM_EvaluatePolynomialInЧебышёвBasisQuantity(benchmark::State& state) {
  std::mt19937_64 random(42);
  std::array<Length, degree + 1> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients[i] = static_cast<double>(random()) * Metre;
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  PolynomialInЧебышёвBasis<Length, Instant, degree> const series(
      coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  Length result = 0.0 * Metre;

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series(t);
      t += Δt;
    }
  }

  // This weird call to `SetLabel` has no effect except that it uses `result`
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

template<int degree>
void BM_EvaluatePolynomialInЧебышёвBasisR3ElementDouble(
    benchmark::State& state) {
  std::mt19937_64 random(42);
  std::array<R3Element<double>, degree + 1> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients[i] = {static_cast<double>(random()),
                       static_cast<double>(random()),
                       static_cast<double>(random())};
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  PolynomialInЧебышёвBasis<R3Element<double>, Instant, degree> const series(
      coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  R3Element<double> result{0.0, 0.0, 0.0};

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series(t);
      t += Δt;
    }
  }

  // This weird call to `SetLabel` has no effect except that it uses `result`
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

template<int degree>
void BM_EvaluatePolynomialInЧебышёвBasisVectorDouble(
    benchmark::State& state) {
  std::mt19937_64 random(42);
  std::array<Multivector<double, ICRS, 1>, degree + 1> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients[i] =
        Multivector<double, ICRS, 1>({static_cast<double>(random()),
                                      static_cast<double>(random()),
                                      static_cast<double>(random())});
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  PolynomialInЧебышёвBasis<Multivector<double, ICRS, 1>, Instant, degree> const
      series(coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  Multivector<double, ICRS, 1> result{};

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series(t);
      t += Δt;
    }
  }

  // This weird call to `SetLabel` has no effect except that it uses `result`
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

template<int degree>
void BM_EvaluatePolynomialInЧебышёвBasisDisplacement(
    benchmark::State& state) {
  std::mt19937_64 random(42);
  std::array<Displacement<ICRS>, degree + 1> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients[i] =
        Displacement<ICRS>({static_cast<double>(random()) * Metre,
                            static_cast<double>(random()) * Metre,
                            static_cast<double>(random()) * Metre});
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  PolynomialInЧебышёвBasis<Displacement<ICRS>, Instant, degree> const series(
      coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  Displacement<ICRS> result{};

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series(t);
      t += Δt;
    }
  }

  // This weird call to `SetLabel` has no effect except that it uses `result`
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDouble, 4)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDouble, 8)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDouble, 15)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDouble, 16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDouble, 17)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDouble, 18)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDouble, 19)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisQuantity, 4)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisQuantity, 8)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisQuantity, 15)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisQuantity, 16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisQuantity, 17)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisQuantity, 18)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisQuantity, 19)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisR3ElementDouble, 4)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisR3ElementDouble, 8)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisR3ElementDouble, 15)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisR3ElementDouble, 16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisR3ElementDouble, 17)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisR3ElementDouble, 18)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisR3ElementDouble, 19)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisVectorDouble, 4)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisVectorDouble, 8)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisVectorDouble, 15)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisVectorDouble, 16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisVectorDouble, 17)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisVectorDouble, 18)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisVectorDouble, 19)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDisplacement, 4)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDisplacement, 8)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDisplacement, 15)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDisplacement, 16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDisplacement, 17)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDisplacement, 18)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInЧебышёвBasisDisplacement, 19)
    ->Unit(benchmark::kMicrosecond);

}  // namespace numerics
}  // namespace principia
