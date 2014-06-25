// .\Release\clr_benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30 --benchmark_filter=HarmonicOscillator  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/16-22:12:38
// Benchmark                               Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------
// BM_CLR_SolveHarmonicOscillator        6683419443 6645642600          1
// BM_CLR_SolveHarmonicOscillator        6677659250 6661242700          1
// BM_CLR_SolveHarmonicOscillator        6656656004 6645642600          1
// BM_CLR_SolveHarmonicOscillator        6658656813 6661242700          1
// BM_CLR_SolveHarmonicOscillator        6639654375 6598842300          1
// BM_CLR_SolveHarmonicOscillator_mean   6663209177 6642522580          1
// BM_CLR_SolveHarmonicOscillator_stddev   15721045   22927371          0

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::clr_benchmarks_adapter::SPRKIntegratorCLRBenchmark;

namespace principia {
namespace clr_benchmarks {

void BM_CLR_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    SPRKIntegratorCLRBenchmark::SolveHarmonicOscillator();
  }
}

BENCHMARK(BM_CLR_SolveHarmonicOscillator);

}  // namespace clr_benchmarks
}  // namespace principia
