
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=Elbdj  // NOLINT(whitespace/line_length)

#include <random>
#include <vector>

#include "benchmark/benchmark.h"
#include "numerics/elliptic_integrals.hpp"
#include "quantities/numbers.hpp"

namespace principia {
namespace numerics {

void BM_Elbdj(benchmark::State& state) {
  constexpr int size = 20;

  std::mt19937_64 random(42);
  // NOTE(phl): A negative lower bound here causes a failure "(elsbdj) too many
  // iterations".  Not sure if this is normal.
  std::uniform_real_distribution<> distribution(0.0, 10.0);
  std::vector<double> phis;
  std::vector<double> ns;
  std::vector<double> mcs;
  for (int i = 0; i < size; ++i) {
    phis.push_back(distribution(random));
    ns.push_back(distribution(random));
    mcs.push_back(distribution(random));
  }

  while (state.KeepRunningBatch(size * size * size)) {
    double b;
    double d;
    double j;
    for (double const phi : phis) {
      double const phic = π / 2 - phi;
      for (double const n : ns) {
        for (double const mc : mcs) {
          Elbdj(phi, phic, n, mc, b, d, j);
        }
      }
    }
    benchmark::DoNotOptimize(b);
    benchmark::DoNotOptimize(d);
    benchmark::DoNotOptimize(j);
  }
}

BENCHMARK(BM_Elbdj);

}  // namespace numerics
}  // namespace principia
