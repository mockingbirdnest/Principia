// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=BM_EvaluateElementaryFunction  // NOLINT(whitespace/line_length)

#include <cmath>
#include <random>

#include "benchmark/benchmark.h"
#include "benchmarks/metric.hpp"
#include "functions/cos.hpp"
#include "functions/sin.hpp"
#include "quantities/numbers.hpp"  // 🧙 For π.

namespace principia {
namespace functions {

using namespace principia::benchmarks::_metric;
using namespace principia::functions::_cos;
using namespace principia::functions::_sin;

static constexpr std::int64_t number_of_iterations = 1000;

template<Metric metric, double (__cdecl *fn)(double)>
void BM_EvaluateElementaryFunction(benchmark::State& state) {
  using Value = double;
  using Argument = double;

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(0.0, π / 4);

  Argument a[number_of_iterations];
  for (std::int64_t i = 0; i < number_of_iterations; ++i) {
    a[i] = uniformly_at(random);
  }

  if constexpr (metric == Metric::Throughput) {
    Value v[number_of_iterations];
    while (state.KeepRunningBatch(number_of_iterations)) {
      for (std::int64_t i = 0; i < number_of_iterations; ++i) {
        v[i] = fn(a[i]);
      }
      benchmark::DoNotOptimize(v);
    }
  } else {
    static_assert(metric == Metric::Latency);
    Value v;
    while (state.KeepRunningBatch(number_of_iterations)) {
      Argument argument = a[0];
      for (std::int64_t i = 0; i < number_of_iterations; ++i) {
        v = fn(argument);
        argument = (v + a[i]) - v;
      }
    }
  }
}

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, std::sin)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, std::sin)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, cr_sin)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, cr_sin)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, std::cos)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, std::cos)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, cr_cos)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, cr_cos)
    ->Unit(benchmark::kNanosecond);

}  // namespace functions
}  // namespace principia
