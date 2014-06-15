
// .\Release\clr_benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/15-10:11:09
// Benchmark                   Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------
// BM_CLR_SolarSystem        33467878796 33150212500          1
// BM_CLR_SolarSystem        34258323646 34086218500          1
// BM_CLR_SolarSystem        32034199215 31933404700          1
// BM_CLR_SolarSystem        31540147997 31387401200          1
// BM_CLR_SolarSystem        34020395256 33930217500          1
// BM_CLR_SolarSystem_mean   33064188982 32897490880          1
// BM_CLR_SolarSystem_stddev 1085065283 1072688104          0
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
