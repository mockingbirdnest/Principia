
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=FastSinCos  // NOLINT(whitespace/line_length)

#include "numerics/fast_sin_cos_cycle.hpp"

#include <random>
#include <vector>

#include "benchmark/benchmark.h"

namespace principia {
namespace numerics {

void BM_FastSinCosCycle(benchmark::State& state) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> const distribution(-1.0, 1.0);
  std::vector<double> input;
  for (int i = 0; i < 1e3; ++i) {
    input.push_back(distribution(random));
  }

  while (state.KeepRunning()) {
    double sin;
    double cos;
    for (double const x : input) {
      FastSinCosCycle(x, sin, cos);
    }
    benchmark::DoNotOptimize(sin);
    benchmark::DoNotOptimize(cos);
  }
}

BENCHMARK(BM_FastSinCosCycle);

}  // namespace numerics
}  // namespace principia
