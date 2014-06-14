
// .\Release\clr_benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/14-19:59:24
// Benchmark                   Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------
// BM_CLR_SolarSystem        31568664766 31481001800          1
// BM_CLR_SolarSystem        32741996742 32697809600          1
// BM_CLR_SolarSystem        32003994919 31933404700          1
// BM_CLR_SolarSystem        31928994633 31871004300          1
// BM_CLR_SolarSystem        32342993237 32307807100          1
// BM_CLR_SolarSystem_mean   32117328859 32058205500          1
// BM_CLR_SolarSystem_stddev  397613135  413564506          0
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
