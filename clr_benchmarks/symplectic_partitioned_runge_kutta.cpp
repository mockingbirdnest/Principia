
// .\Release\clr_benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30 --benchmark_filter=HarmonicOscillator  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/09-00:26:58
// Benchmark                               Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------
// BM_CLR_SolveHarmonicOscillator        3304006077 3291621100          2
// BM_CLR_SolveHarmonicOscillator        3268822528 3260420900          2
// BM_CLR_SolveHarmonicOscillator        3269823409 3268220950          2
// BM_CLR_SolveHarmonicOscillator        3269323803 3260420900          2
// BM_CLR_SolveHarmonicOscillator        3272821642 3276021000          2
// BM_CLR_SolveHarmonicOscillator_mean   3276959492 3271340970          2
// BM_CLR_SolveHarmonicOscillator_stddev   13594684   11674046          0

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
