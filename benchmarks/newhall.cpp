
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=Newhall  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/05/24-13:16:32
// Benchmark                    Time(ns)    CPU(ns) Iterations
// -----------------------------------------------------------
// BM_NewhallApproximation/4         589        562    2000000
// BM_NewhallApproximation/8         657        624    2000000
// BM_NewhallApproximation/16        754        741    2000000

#include <random>
#include <vector>

#include "benchmark/benchmark.h"
#include "geometry/named_quantities.hpp"
#include "numerics/newhall.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::Instant;
using quantities::Variation;
using quantities::si::Second;

namespace numerics {

void BM_NewhallApproximation(benchmark::State& state) {
  int const degree = state.range_x();
  std::mt19937_64 random(42);
  std::vector<double> p;
  std::vector<Variation<double>> v;
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;

  while (state.KeepRunning()) {
    state.PauseTiming();
    p.clear();
    v.clear();
    for (int i = 0; i <= 8; ++i) {
      p.push_back(static_cast<double>(static_cast<double>(random())));
      v.push_back(static_cast<double>(static_cast<double>(random())) / Second);
    }
    state.ResumeTiming();
    auto const series =
        NewhallApproximationInЧебышёвBasis<double>(degree, p, v, t_min, t_max);
  }
}

BENCHMARK(BM_NewhallApproximation)->Arg(4)->Arg(8)->Arg(16);

}  // namespace numerics
}  // namespace principia
