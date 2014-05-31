#define GLOG_NO_ABBREVIATED_SEVERITIES

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

// using principia::clr_benchmarks::SPRKIntegratorCLRBenchmark;

namespace principia {
namespace benchmarks {

static void BM_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    principia::clr_benchmarks::SPRKIntegratorCLRBenchmark::SolveHarmonicOscillator();
  }
}

BENCHMARK(BM_SolveHarmonicOscillator);

}  // namespace benchmarks
}  // namespace principia
