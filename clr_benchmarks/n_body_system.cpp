
// .\Release\clr_benchmarks.exe --benchmark_filter=Solar --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/10-05:17:42
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_CLR_SimulateSolarSystem        49901121345 49639518200          1
// BM_CLR_SimulateSolarSystem        47652830568 47377503700          1
// BM_CLR_SimulateSolarSystem        49001193014 48765912600          1
// BM_CLR_SimulateSolarSystem        54341507879 54054346500          1
// BM_CLR_SimulateSolarSystem        48524459091 47954707400          1
// BM_CLR_SimulateSolarSystem_mean   49884222379 49558397680          5
// BM_CLR_SimulateSolarSystem_stddev 2344162370 2373869836          5

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
