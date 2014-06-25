// .\Release\clr_benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/15-23:53:12
// Benchmark                   Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------
// BM_CLR_SolarSystem        32835119701 32791410200          1
// BM_CLR_SolarSystem        32463242746 32339007300          1
// BM_CLR_SolarSystem        32606255270 32557408700          1
// BM_CLR_SolarSystem        35357531525 35053424700          1
// BM_CLR_SolarSystem        35877581812 35802229500          1
// BM_CLR_SolarSystem_mean   33827946211 33708696080          1
// BM_CLR_SolarSystem_stddev 1475214297 1430671214          0
#include "benchmark/benchmark.h"

using principia::clr_benchmarks_adapter::NBodySystemCLRBenchmark;

namespace principia {
namespace clr_benchmarks {

void BM_CLR_SolarSystem(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    NBodySystemCLRBenchmark::SimulateSolarSystem();
  }
}

BENCHMARK(BM_CLR_SolarSystem);

}  // namespace clr_benchmarks
}  // namespace principia
