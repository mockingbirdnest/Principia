
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=Gscd  // NOLINT(whitespace/line_length)

#include <random>
#include <vector>

#include "benchmark/benchmark.h"
#include "numerics/elliptic_functions.hpp"
#include "quantities/numbers.hpp"

namespace principia {
namespace numerics {

void BM_Gscd(benchmark::State& state) {
  constexpr int size = 100;

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> distribution_u(-10.0, 10.0);
  std::uniform_real_distribution<> distribution_mc(0.0, 1.0);
  std::vector<double> us;
  std::vector<double> mcs;
  for (int i = 0; i < size; ++i) {
    us.push_back(distribution_u(random));
    mcs.push_back(distribution_mc(random));
  }

  while (state.KeepRunningBatch(size * size)) {
    double s;
    double c;
    double d;
    for (double const u : us) {
      for (double const mc : mcs) {
        Gscd(u, mc, s, c, d);
      }
    }
    benchmark::DoNotOptimize(s);
    benchmark::DoNotOptimize(c);
    benchmark::DoNotOptimize(d);
  }
}

BENCHMARK(BM_Gscd);

}  // namespace numerics
}  // namespace principia