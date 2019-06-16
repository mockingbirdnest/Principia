
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=FukushimaEllipticBDJ  // NOLINT(whitespace/line_length)

#include <random>
#include <vector>

#include "benchmark/benchmark.h"
#include "numerics/elliptic_integrals.hpp"
#include "quantities/numbers.hpp"

namespace principia {
namespace numerics {

void BM_FukushimaEllipticBDJ(benchmark::State& state) {
  constexpr int size = 20;

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> distribution_φ(0.0, π / 2);
  std::uniform_real_distribution<> distribution_n(0.0, 1.0);
  std::uniform_real_distribution<> distribution_mc(0.0, 1.0);
  std::vector<double> φs;
  std::vector<double> ns;
  std::vector<double> mcs;
  for (int i = 0; i < size; ++i) {
    φs.push_back(distribution_φ(random));
    ns.push_back(distribution_n(random));
    mcs.push_back(distribution_mc(random));
  }

  while (state.KeepRunningBatch(size * size * size)) {
    double b;
    double d;
    double j;
    for (double const φ : φs) {
      for (double const n : ns) {
        for (double const mc : mcs) {
          FukushimaEllipticBDJ(φ, n, mc, b, d, j);
        }
      }
    }
    benchmark::DoNotOptimize(b);
    benchmark::DoNotOptimize(d);
    benchmark::DoNotOptimize(j);
  }
}

BENCHMARK(BM_FukushimaEllipticBDJ);

}  // namespace numerics
}  // namespace principia
