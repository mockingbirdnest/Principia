
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=Elliptic  // NOLINT(whitespace/line_length)

#include <random>
#include <vector>

#include "benchmark/benchmark.h"
#include "numerics/elliptic_integrals.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"

namespace principia {

using quantities::Angle;
using quantities::si::Radian;

namespace numerics {

void BM_EllipticEFΠ(benchmark::State& state) {
  constexpr int size = 20;

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> distribution_φ(0.0, π / 2);
  std::uniform_real_distribution<> distribution_n(0.0, 1.0);
  std::uniform_real_distribution<> distribution_mc(0.0, 1.0);
  std::vector<Angle> φs;
  std::vector<double> ns;
  std::vector<double> mcs;
  for (int i = 0; i < size; ++i) {
    φs.push_back(distribution_φ(random) * Radian);
    ns.push_back(distribution_n(random));
    mcs.push_back(distribution_mc(random));
  }

  while (state.KeepRunningBatch(size * size * size)) {
    double e;
    double f;
    double ᴨ;
    for (Angle const φ : φs) {
      for (double const n : ns) {
        for (double const mc : mcs) {
          EllipticEFΠ(φ, n, mc, e, f, ᴨ);
        }
      }
    }
    benchmark::DoNotOptimize(e);
    benchmark::DoNotOptimize(f);
    benchmark::DoNotOptimize(ᴨ);
  }
}

void BM_FukushimaEllipticBDJ(benchmark::State& state) {
  constexpr int size = 20;

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> distribution_φ(0.0, π / 2);
  std::uniform_real_distribution<> distribution_n(0.0, 1.0);
  std::uniform_real_distribution<> distribution_mc(0.0, 1.0);
  std::vector<Angle> φs;
  std::vector<double> ns;
  std::vector<double> mcs;
  for (int i = 0; i < size; ++i) {
    φs.push_back(distribution_φ(random) * Radian);
    ns.push_back(distribution_n(random));
    mcs.push_back(distribution_mc(random));
  }

  while (state.KeepRunningBatch(size * size * size)) {
    double b;
    double d;
    double j;
    for (Angle const φ : φs) {
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

BENCHMARK(BM_EllipticEFΠ);
BENCHMARK(BM_FukushimaEllipticBDJ);

}  // namespace numerics
}  // namespace principia
