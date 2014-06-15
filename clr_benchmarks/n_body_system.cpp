
// .\Release\clr_benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/15-15:51:05
// Benchmark                   Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------
// BM_CLR_SolarSystem        31901581264 31808603900          1
// BM_CLR_SolarSystem        32663262583 32557408700          1
// BM_CLR_SolarSystem        31907186485 31808603900          1
// BM_CLR_SolarSystem        31813779808 31793003800          1
// BM_CLR_SolarSystem        34345427967 34195419200          1
// BM_CLR_SolarSystem_mean   32526247621 32432607900          1
// BM_CLR_SolarSystem_stddev  960129716  928540178          0
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
