
//  .\Release\clr_benchmarks.exe  --benchmark_filter=Solar --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/11-21:50:34
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_CLR_SimulateSolarSystem        31866679646 31855404200          1
// BM_CLR_SimulateSolarSystem        33116308190 32978611400          1
// BM_CLR_SimulateSolarSystem        33756355263 33290613400          1
// BM_CLR_SimulateSolarSystem        32331227518 32245406700          1
// BM_CLR_SimulateSolarSystem        33322325445 33321813600          1
// BM_CLR_SimulateSolarSystem_mean   32878579212 32738369860          1
// BM_CLR_SimulateSolarSystem_stddev  685193577  587492393          0
#include "benchmark/benchmark.h"

using principia::clr_benchmarks_adapter::NBodySystemCLRBenchmark;

namespace principia {
namespace clr_benchmarks {

void BM_CLR_SimulateSolarSystem(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    NBodySystemCLRBenchmark::SimulateSolarSystem();
  }
}

BENCHMARK(BM_CLR_SimulateSolarSystem);

}  // namespace clr_benchmarks
}  // namespace principia
