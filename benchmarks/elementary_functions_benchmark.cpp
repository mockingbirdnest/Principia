// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=BM_EvaluateElementaryFunction  // NOLINT(whitespace/line_length)

#include <cmath>
#include <random>

#include "benchmark/benchmark.h"

namespace principia {
namespace functions {

static constexpr std::int64_t number_of_iterations = 1000;

enum class Metric { Latency, Throughput };

template<Metric metric, double (__cdecl *fn)(double)>
void BM_EvaluateElementaryFunction(benchmark::State& state) {
  using Value = double;
  using Argument = double;

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(-10.0, 10.0);
  Argument argument = uniformly_at(random);

  if constexpr (metric == Metric::Throughput) {
    Value v[number_of_iterations];
    while (state.KeepRunningBatch(number_of_iterations)) {
      for (std::int64_t i = 0; i < number_of_iterations; ++i) {
        v[i] = fn(argument);
      }
      benchmark::DoNotOptimize(v);
    }
  } else {
    static_assert(metric == Metric::Latency);
    Value v;
    while (state.KeepRunningBatch(number_of_iterations)) {
      for (std::int64_t i = 0; i < number_of_iterations; ++i) {
        v = fn(argument);
        argument = v;
      }
    }
  }
}

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, std::sin)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, std::sin)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, std::cos)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, std::cos)
    ->Unit(benchmark::kNanosecond);

}  // namespace functions
}  // namespace principia
