
// .\Release\clr_benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/01-20:19:58
// Benchmark                               Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------
// BM_CLR_SolveHarmonicOscillator        2877869106 2880818467          3
// BM_CLR_SolveHarmonicOscillator        2868951775 2865218367          3
// BM_CLR_SolveHarmonicOscillator        2871617294 2875618433          3
// BM_CLR_SolveHarmonicOscillator        2872283515 2870418400          3
// BM_CLR_SolveHarmonicOscillator        2864615815 2865218367          3
// BM_CLR_SolveHarmonicOscillator_mean   2871067501 2871458407         15
// BM_CLR_SolveHarmonicOscillator_stddev    4339202    6064229         15

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
