
// .\Release\clr_benchmarks.exe --benchmark_filter=Solar --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/10-04:19:13
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_CLR_SimulateSolarSystem        42617064443 42385471700          1
// BM_CLR_SimulateSolarSystem        42213265863 41823868100          1
// BM_CLR_SimulateSolarSystem        39901499239 39717854600          1
// BM_CLR_SimulateSolarSystem        43798215595 43305877600          1
// BM_CLR_SimulateSolarSystem        41949605194 41683467200          1
// BM_CLR_SimulateSolarSystem_mean   42095930067 41783307840          5
// BM_CLR_SimulateSolarSystem_stddev 1266496141 1179705926          5

#include "benchmark/benchmark.h"

using principia::clr_benchmarks_adapter::NBodySystemCLRBenchmark;

namespace principia {
namespace clr_benchmarks {

static void BM_CLR_SimulateSolarSystem(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    NBodySystemCLRBenchmark::SimulateSolarSystem();
  }
}

BENCHMARK(BM_CLR_SimulateSolarSystem);

}  // namespace clr_benchmarks
}  // namespace principia
