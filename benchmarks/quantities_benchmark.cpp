// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=DiscreteCosine  // NOLINT(whitespace/line_length)

#include<vector>

#include "benchmark/benchmark.h"
#include "benchmarks/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace quantities {

using namespace principia::benchmarks::_quantities;
using namespace principia::quantities::_named_quantities;

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
