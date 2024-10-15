// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=BM_EvaluateElementaryFunction  // NOLINT(whitespace/line_length)

#include <cmath>
#include <random>

#include "absl/strings/str_cat.h"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_REPEAT.
#include "benchmark/benchmark.h"
#include "benchmarks/metric.hpp"
#include "core-math/cos.h"
#include "core-math/sin.h"
#include "numerics/sin_cos.hpp"
#include "quantities/numbers.hpp"  // 🧙 For π.

namespace principia {
namespace functions {

using namespace principia::benchmarks::_metric;
using namespace principia::numerics::_sin_cos;

static constexpr std::int64_t number_of_iterations = 1000;

template<Metric metric, double (__cdecl *fn)(double)>
void BM_EvaluateElementaryFunction(benchmark::State& state) {
  using Value = double;
  using Argument = double;

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(-2 * π, 2 * π);

  Argument a[number_of_iterations];
  for (std::int64_t i = 0; i < number_of_iterations; ++i) {
    a[i] = uniformly_at(random);
  }

  std::uint64_t iterations = 0;
  std::uint64_t cycles = 0;
  if constexpr (metric == Metric::Throughput) {
    Value v[number_of_iterations];
    while (state.KeepRunningBatch(number_of_iterations)) {
      std::uint64_t const start = __rdtsc();
      for (std::int64_t i = 0; i < number_of_iterations;) {
        PRINCIPIA_REPEAT8(v[i] = fn(a[i]); ++i;)
      }
      std::uint64_t const stop = __rdtsc();
      ++iterations;
      cycles += stop - start;
    }
    benchmark::DoNotOptimize(v);
  } else {
    static_assert(metric == Metric::Latency);
    Value v;
    while (state.KeepRunningBatch(number_of_iterations)) {
      std::uint64_t const start = __rdtsc();
      Argument argument = a[number_of_iterations - 1];
      for (std::int64_t i = 0; i < number_of_iterations; ++i) {
        v = fn(argument);
        argument = (v + a[i]) - v;
      }
      std::uint64_t const stop = __rdtsc();
      ++iterations;
      cycles += stop - start;
    }
    benchmark::DoNotOptimize(v);
  }
  state.SetLabel(
      absl::StrCat("cycles: ",
                   static_cast<double>(cycles) /
                       static_cast<double>(iterations * number_of_iterations)));
}

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, std::sin)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, std::sin)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, cr_sin)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, cr_sin)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, Sin)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, Sin)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, std::cos)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, std::cos)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, cr_cos)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, cr_cos)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Latency, Cos)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateElementaryFunction, Metric::Throughput, Cos)
    ->Unit(benchmark::kNanosecond);

}  // namespace functions
}  // namespace principia
