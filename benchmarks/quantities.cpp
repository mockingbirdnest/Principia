
#include "benchmarks/quantities.hpp"

#include<vector>

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::base::check_not_null;

namespace principia {
namespace benchmarks {

void BM_DimensionfulDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)>
  std::vector<quantities::Momentum> output;
  while (state.KeepRunning()) {
    DimensionfulDiscreteCosineTransform(&output);
  }
}
BENCHMARK(BM_DimensionfulDiscreteCosineTransform);

void BM_DoubleDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)>
  std::vector<double> output;
  while (state.KeepRunning()) {
    DoubleDiscreteCosineTransform(&output);
  }
}
BENCHMARK(BM_DoubleDiscreteCosineTransform);

}  // namespace benchmarks
}  // namespace principia
