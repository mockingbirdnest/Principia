
// .\Release\clr_benchmarks.exe  --benchmark_filter=Solar --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/11-22:02:54
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_CLR_SolarSystem        31282008148 31278200500          1
// BM_CLR_SolarSystem        32160996176 32136206000          1
// BM_CLR_SolarSystem        32204995848 32198606400          1
// BM_CLR_SolarSystem        32454995132 32448208000          1
// BM_CLR_SolarSystem        31222993589 31168999800          1
// BM_CLR_SolarSystem_mean   31865197779 31846044140          1
// BM_CLR_SolarSystem_stddev  510559831  519986693          0
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
