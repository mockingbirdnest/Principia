
#include "benchmarks/quantities.hpp"

#include<vector>

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {
namespace benchmarks {

static void BM_DimensionfulDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)>
  std::vector<quantities::Momentum> output;
  while (state.KeepRunning()) {
    DimensionfulDiscreteCosineTransform(&output);
  }
}
BENCHMARK(BM_DimensionfulDiscreteCosineTransform);

static void BM_DoubleDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)>
  std::vector<double> output;
  while (state.KeepRunning()) {
    DoubleDiscreteCosineTransform(&output);
  }
}
BENCHMARK(BM_DoubleDiscreteCosineTransform);

}  // namespace benchmarks
}  // namespace principia
