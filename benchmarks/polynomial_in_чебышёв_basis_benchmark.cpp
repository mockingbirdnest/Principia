// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Evaluate  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/06/09-22:18:17
// Benchmark                    Time(ns)    CPU(ns) Iterations
// -----------------------------------------------------------
// BM_EvaluateDouble/4              8031       8139      63254
// BM_EvaluateDouble/8             15548      15602      32995
// BM_EvaluateDouble/15            35945      36121      14252
// BM_EvaluateDouble/16            36019      36104      14259
// BM_EvaluateDouble/17            39483      39558      13014
// BM_EvaluateDouble/18            42357      42413      12138
// BM_EvaluateDouble/19            45989      44785      11495
// BM_EvaluateDisplacement/4       52921      53034       9707
// BM_EvaluateDisplacement/8      101641     101800       5057
// BM_EvaluateDisplacement/15     186835     187269       2749
// BM_EvaluateDisplacement/16     199863     200781       2564
// BM_EvaluateDisplacement/17     211855     212203       2426
// BM_EvaluateDisplacement/18     225181     225494       2283
// BM_EvaluateDisplacement/19     236922     237237       2170

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
void BM_EvaluateDouble(benchmark::State& state) {
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

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  state.SetLabel(std::to_string(result).substr(0, 0));
}

template<int degree>
void BM_EvaluateQuantity(benchmark::State& state) {
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

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

template<int degree>
void BM_EvaluateR3ElementDouble(benchmark::State& state) {
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

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

template<int degree>
void BM_EvaluateVectorDouble(benchmark::State& state) {
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

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

template<int degree>
void BM_EvaluateDisplacement(benchmark::State& state) {
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

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

BENCHMARK_TEMPLATE(BM_EvaluateDouble, 4)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDouble, 8)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDouble, 15)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDouble, 16)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDouble, 17)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDouble, 18)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDouble, 19)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_EvaluateQuantity, 4)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateQuantity, 8)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateQuantity, 15)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateQuantity, 16)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateQuantity, 17)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateQuantity, 18)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateQuantity, 19)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_EvaluateR3ElementDouble, 4)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateR3ElementDouble, 8)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateR3ElementDouble, 15)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateR3ElementDouble, 16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateR3ElementDouble, 17)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateR3ElementDouble, 18)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateR3ElementDouble, 19)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_EvaluateVectorDouble, 4)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateVectorDouble, 8)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateVectorDouble, 15)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateVectorDouble, 16)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateVectorDouble, 17)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateVectorDouble, 18)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateVectorDouble, 19)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_EvaluateDisplacement, 4)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDisplacement, 8)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDisplacement, 15)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDisplacement, 16)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDisplacement, 17)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDisplacement, 18)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EvaluateDisplacement, 19)->Unit(benchmark::kMicrosecond);

}  // namespace numerics
}  // namespace principia
