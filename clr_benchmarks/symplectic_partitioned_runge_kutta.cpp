#define GLOG_NO_ABBREVIATED_SEVERITIES

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {
namespace clr_benchmarks {

static void BM_CLR_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    principia::clr_benchmarks_adapter::SPRKIntegratorCLRBenchmark::SolveHarmonicOscillator();
  }
}

BENCHMARK(BM_CLR_SolveHarmonicOscillator);

}  // namespace clr_benchmarks
}  // namespace principia
