
#include "benchmarks/quantities.hpp"

#include<vector>

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {
namespace quantities {

void BM_DimensionfulDiscreteCosineTransform(benchmark::State& state) {
  std::vector<Momentum> output;
  while (state.KeepRunning()) {
    DimensionfulDiscreteCosineTransform(output);
  }
}
BENCHMARK(BM_DimensionfulDiscreteCosineTransform);

void BM_DoubleDiscreteCosineTransform(benchmark::State& state) {
  std::vector<double> output;
  while (state.KeepRunning()) {
    DoubleDiscreteCosineTransform(output);
  }
}
BENCHMARK(BM_DoubleDiscreteCosineTransform);

}  // namespace quantities
}  // namespace principia
