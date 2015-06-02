
// .\Release\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Evaluate  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/05/15-18:21:37
// Benchmark              Time(ns)    CPU(ns) Iterations
// -----------------------------------------------------
// BM_EvaluateDouble/4       11859      12009      42869
// BM_EvaluateDouble/8       18978      19040      27038
// BM_EvaluateDouble/15      40024      40119      12832
// BM_EvaluateDouble/16      38326      38407      13404
// BM_EvaluateDouble/17      43830      44042      11689
// BM_EvaluateDouble/18      49399      49519      10396
// BM_EvaluateDouble/19      54013      54121       9512

// .\Release\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=Newhall  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/05/24-13:16:32
// Benchmark                    Time(ns)    CPU(ns) Iterations
// -----------------------------------------------------------
// BM_NewhallApproximation/4         589        562    2000000
// BM_NewhallApproximation/8         657        624    2000000
// BM_NewhallApproximation/16        754        741    2000000

#include <random>
#include <vector>

#include "quantities/si.hpp"
#include "numerics/чебышёв_series.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using numerics::ЧебышёвSeries;
using si::Second;

namespace benchmarks {

namespace {
int const kEvaluationsPerIteration = 1000;
}  // namespace

void BM_EvaluateDouble(benchmark::State& state) {  // NOLINT(runtime/references)
  state.PauseTiming();
  int const degree = state.range_x();
  std::mt19937_64 random(42);
  std::vector<double> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients.push_back(static_cast<double>(random()));
  }
  Instant const t_min(random() * Second);
  Instant const t_max = t_min + random() * Second;
  ЧебышёвSeries<double> const series(coefficients, t_min, t_max);

  Instant t = t_min;
  Time const ∆t = (t_max - t_min) * 1E-9;
  double result = 0.0;

  state.ResumeTiming();
  while (state.KeepRunning()) {
    for (int i = 0; i < kEvaluationsPerIteration; ++i) {
      result += series.Evaluate(t);
      t += ∆t;
    }
  }
  state.PauseTiming();

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  state.SetLabel(std::to_string(result).substr(0, 0));
}

void BM_NewhallApproximation(
    benchmark::State& state) {  // NOLINT(runtime/references)
  state.PauseTiming();
  int const degree = state.range_x();
  std::mt19937_64 random(42);
  std::vector<double> p;
  std::vector<Variation<double>> v;
  Instant const t_min(random() * Second);
  Instant const t_max = t_min + random() * Second;

  while (state.KeepRunning()) {
    p.clear();
    v.clear();
    for (int i = 0; i <= 8; ++i) {
      p.push_back(static_cast<double>(random()));
      v.push_back(static_cast<double>(random()) / Second);
    }
    state.ResumeTiming();
    auto const series =
        ЧебышёвSeries<double>::NewhallApproximation(degree, p, v, t_min, t_max);
    state.PauseTiming();
  }
}

BENCHMARK(BM_EvaluateDouble)->
    Arg(4)->Arg(8)->Arg(15)->Arg(16)->Arg(17)->Arg(18)->Arg(19);
BENCHMARK(BM_NewhallApproximation)->
    Arg(4)->Arg(8)->Arg(16);

}  // namespace benchmarks
}  // namespace principia
