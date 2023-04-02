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
#include <vector>

#include "astronomy/frames.hpp"
#include "benchmark/benchmark.h"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "numerics/чебышёв_series.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {

using namespace principia::astronomy::_frames;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_space;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

namespace {
constexpr int evaluations_per_iteration = 1000;
}  // namespace

void BM_EvaluateDouble(benchmark::State& state) {
  int const degree = state.range(0);
  std::mt19937_64 random(42);
  std::vector<double> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients.push_back(static_cast<double>(random()));
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  ЧебышёвSeries<double> const series(coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  double result = 0.0;

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series.Evaluate(t);
      t += Δt;
    }
  }

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  state.SetLabel(std::to_string(result).substr(0, 0));
}

void BM_EvaluateQuantity(benchmark::State& state) {
  int const degree = state.range(0);
  std::mt19937_64 random(42);
  std::vector<Length> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients.push_back(static_cast<double>(random()) * Metre);
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  ЧебышёвSeries<Length> const series(coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  Length result = 0.0 * Metre;

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series.Evaluate(t);
      t += Δt;
    }
  }

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

void BM_EvaluateR3ElementDouble(benchmark::State& state) {
  int const degree = state.range(0);
  std::mt19937_64 random(42);
  std::vector<R3Element<double>> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients.push_back({static_cast<double>(random()),
                            static_cast<double>(random()),
                            static_cast<double>(random())});
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  ЧебышёвSeries<R3Element<double>> const series(coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  R3Element<double> result{0.0, 0.0, 0.0};

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series.Evaluate(t);
      t += Δt;
    }
  }

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

void BM_EvaluateVectorDouble(benchmark::State& state) {
  int const degree = state.range(0);
  std::mt19937_64 random(42);
  std::vector<Multivector<double, ICRS, 1>> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients.push_back(
        Multivector<double, ICRS, 1>({static_cast<double>(random()),
                                      static_cast<double>(random()),
                                      static_cast<double>(random())}));
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  ЧебышёвSeries<Multivector<double, ICRS, 1>> const series(
      coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  Multivector<double, ICRS, 1> result{};

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series.Evaluate(t);
      t += Δt;
    }
  }

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

void BM_EvaluateDisplacement(benchmark::State& state) {
  int const degree = state.range(0);
  std::mt19937_64 random(42);
  std::vector<Displacement<ICRS>> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients.push_back(
        Displacement<ICRS>({static_cast<double>(random()) * Metre,
                            static_cast<double>(random()) * Metre,
                            static_cast<double>(random()) * Metre}));
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  ЧебышёвSeries<Displacement<ICRS>> const series(coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1e-9;
  Displacement<ICRS> result{};

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += series.Evaluate(t);
      t += Δt;
    }
  }

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

BENCHMARK(BM_EvaluateDouble)
    ->Arg(4)
    ->Arg(8)
    ->Arg(15)
    ->Arg(16)
    ->Arg(17)
    ->Arg(18)
    ->Arg(19)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_EvaluateQuantity)
    ->Arg(4)
    ->Arg(8)
    ->Arg(15)
    ->Arg(16)
    ->Arg(17)
    ->Arg(18)
    ->Arg(19)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_EvaluateR3ElementDouble)
    ->Arg(4)
    ->Arg(8)
    ->Arg(15)
    ->Arg(16)
    ->Arg(17)
    ->Arg(18)
    ->Arg(19)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_EvaluateVectorDouble)
    ->Arg(4)
    ->Arg(8)
    ->Arg(15)
    ->Arg(16)
    ->Arg(17)
    ->Arg(18)
    ->Arg(19)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_EvaluateDisplacement)
    ->Arg(4)
    ->Arg(8)
    ->Arg(15)
    ->Arg(16)
    ->Arg(17)
    ->Arg(18)
    ->Arg(19)
    ->Unit(benchmark::kMicrosecond);

}  // namespace numerics
}  // namespace principia
