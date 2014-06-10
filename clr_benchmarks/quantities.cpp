#include "benchmark/benchmark.h"

using principia::clr_benchmarks_adapter::QuantitiesCLRBenchmark;

namespace principia {
namespace clr_benchmarks {

void BM_CLR_DimensionfulDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    QuantitiesCLRBenchmark::DimensionfulDiscreteCosineTransform();
  }
}

BENCHMARK(BM_CLR_DimensionfulDiscreteCosineTransform);

void BM_CLR_DoubleDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    QuantitiesCLRBenchmark::DoubleDiscreteCosineTransform();
  }
}

BENCHMARK(BM_CLR_DoubleDiscreteCosineTransform);

}  // namespace clr_benchmarks
}  // namespace principia
