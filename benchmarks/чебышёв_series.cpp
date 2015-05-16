// C:\Users\phl\Projects\GitHub\Principia [Series +0 ~1 -0]> .\Release\benchmarks.exe --benchmark_filter=Evaluate  // NOLINT(whitespace/line_length)
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

#include <random>

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

void BM_EvaluateDouble(benchmark::State& state) {
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

BENCHMARK(BM_EvaluateDouble)->
    Arg(4)->Arg(8)->Arg(15)->Arg(16)->Arg(17)->Arg(18)->Arg(19);

}  // namespace benchmarks
}  // namespace principia
