#include "benchmarks/quantities.hpp"

#include<vector>

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {
namespace quantities {

void BM_DimensionfulDiscreteCosineTransform(benchmark::State& state) {
  std::vector<Momentum> output;
  for (auto _ : state) {
    DimensionfulDiscreteCosineTransform(output);
  }
}
BENCHMARK(BM_DimensionfulDiscreteCosineTransform);

void BM_DoubleDiscreteCosineTransform(benchmark::State& state) {
  std::vector<double> output;
  for (auto _ : state) {
    DoubleDiscreteCosineTransform(output);
  }
}
BENCHMARK(BM_DoubleDiscreteCosineTransform);

}  // namespace quantities
}  // namespace principia
