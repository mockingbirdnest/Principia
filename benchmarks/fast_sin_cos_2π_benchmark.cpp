
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=FastSinCos  // NOLINT(whitespace/line_length)

#include "numerics/fast_sin_cos_2π.hpp"

#include <random>
#include <vector>

#include "benchmark/benchmark.h"
#include "quantities/numbers.hpp"

namespace principia {
namespace numerics {
namespace {

// Returns a value which has a data dependency on both cos_2πx and sin_2πx.
// The result is (cos_2πx bitand cos_mask) bitxor sin_2πx.
// The latency is between 1 and 2 cycles: at worst this is an and and an xor, at
// best the xor can only be computed given both trigonometric lines.
double MixTrigonometricLines(double cos_2πx, double sin_2πx, __m128d const cos_mask) {
  __m128d const cos_mantissa_bits =
      _mm_and_pd(_mm_set_sd(cos_2πx), cos_mask);
  __m128d const sin_all_bits = _mm_set_sd(sin_2πx);
  __m128d const mixed_bits = _mm_xor_pd(cos_mantissa_bits, sin_all_bits);
  return _mm_cvtsd_f64(mixed_bits);
}

static const __m128d mantissa_bits =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x000F'FFFF'FFFF'FFFF));
// When iterated, the quadrant of the result is unbiased.
double ThoroughlyMixTrigonometricLines(double cos_2πx, double sin_2πx) {
  return MixTrigonometricLines(cos_2πx, sin_2πx, mantissa_bits);
}

static const __m128d mantissa_bits_and_5_bits_of_exponent =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x01FF'FFFF'FFFF'FFFF));
// Same as above, but when iterated, the result is quickly confined to [0, 1/8].
double PoorlyMixTrigonometricLines(double cos_2πx, double sin_2πx) {
  return MixTrigonometricLines(
      cos_2πx, sin_2πx, mantissa_bits_and_5_bits_of_exponent);
}

}  // namespace

void BM_FastSinCos2πPoorlyPredictedLatency(benchmark::State& state) {
  while (state.KeepRunning()) {
    double sin = π;
    double cos = 0.0;
    for (int i = 0; i < 1e3; ++i) {
      FastSinCos2π(ThoroughlyMixTrigonometricLines(cos, sin), sin, cos);
    }
  }
}

void BM_FastSinCos2πWellPredictedLatency(benchmark::State& state) {
  while (state.KeepRunning()) {
    double sin = π;
    double cos = 0.0;
    for (int i = 0; i < 1e3; ++i) {
      FastSinCos2π(PoorlyMixTrigonometricLines(cos, sin), sin, cos);
    }
  }
}

void BM_FastSinCos2πThroughput(benchmark::State& state) {
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
      FastSinCos2π(x, sin, cos);
    }
    benchmark::DoNotOptimize(sin);
    benchmark::DoNotOptimize(cos);
  }
}

BENCHMARK(BM_FastSinCos2πPoorlyPredictedLatency);
BENCHMARK(BM_FastSinCos2πWellPredictedLatency);
BENCHMARK(BM_FastSinCos2πThroughput);

}  // namespace numerics
}  // namespace principia
