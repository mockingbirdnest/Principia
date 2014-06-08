
// .\Release\clr_benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30 --benchmark_filter=HarmonicOscillator  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/09-00:26:58
// Benchmark                               Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------
// BM_CLR_SolveHarmonicOscillator        3326480366 3330621350          2
// BM_CLR_SolveHarmonicOscillator        3314327217 3315021250          2
// BM_CLR_SolveHarmonicOscillator        3336829401 3338421400          2
// BM_CLR_SolveHarmonicOscillator        3315327741 3315021250          2
// BM_CLR_SolveHarmonicOscillator        3302824355 3299421150          2
// BM_CLR_SolveHarmonicOscillator_mean   3319157816 3319701280          2
// BM_CLR_SolveHarmonicOscillator_stddev   11581803   13599852          0

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::clr_benchmarks_adapter::SPRKIntegratorCLRBenchmark;

namespace principia {
namespace clr_benchmarks {

static void BM_CLR_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    SPRKIntegratorCLRBenchmark::SolveHarmonicOscillator();
  }
}

BENCHMARK(BM_CLR_SolveHarmonicOscillator);

}  // namespace clr_benchmarks
}  // namespace principia
