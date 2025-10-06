#include <emmintrin.h>

#include "nanobenchmarks/function_registry.hpp"  // 🧙 For BENCHMARK_FUNCTION etc.
#include "numerics/cbrt.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _examples {

using namespace principia::numerics::_cbrt;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;

BENCHMARKED_FUNCTION(twice, double, double) {
  return 2 * x;
}

BENCHMARKED_FUNCTION(thrice) {
  return 3 * x;
}

BENCHMARKED_FUNCTION(inc) {
  return x + 1;
}

BENCHMARKED_FUNCTION(multiply_4_times) {
  return x * x * x * x * x;
}

BENCHMARKED_FUNCTION(add_16_times) {
  return x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x;
}

BENCHMARKED_FUNCTION(square_root) {
  __m128d x_0 = _mm_set_sd(x);
  return _mm_cvtsd_f64(_mm_sqrt_sd(x_0, x_0));
}

BENCHMARKED_FUNCTION(sqrt_sqrt) {
  __m128d x_0 = _mm_set_sd(x);
  x_0 = _mm_sqrt_sd(x_0, x_0);
  return _mm_cvtsd_f64(_mm_sqrt_sd(x_0, x_0));
}

BENCHMARKED_FUNCTION(square_root_division) {
  __m128d x_0 = _mm_set_sd(x);
  return _mm_cvtsd_f64(_mm_div_sd(x_0, _mm_sqrt_sd(x_0, x_0)));
}
BENCHMARK_FUNCTION(Cbrt);

using namespace principia::numerics::_cbrt::internal;

BENCHMARK_FUNCTION(method_3²ᴄZ5¹::Cbrt<Rounding::Faithful>);
BENCHMARK_FUNCTION(method_3²ᴄZ5¹::Cbrt<Rounding::Correct>);
#if PRINCIPIA_COMPILER_MSVC
BENCHMARK_FUNCTION(method_5²Z4¹FMA::Cbrt<Rounding::Faithful>);
BENCHMARK_FUNCTION(method_5²Z4¹FMA::Cbrt<Rounding::Correct>);
#endif

BENCHMARKED_FUNCTION(std_sin) {
  return std::sin(x);
}

BENCHMARKED_FUNCTION(principia_sin) {
  return Sin(x * Radian);
}

BENCHMARKED_FUNCTION(std_cos) {
  return std::cos(x);
}

BENCHMARKED_FUNCTION(principia_cos) {
  return Cos(x * Radian);
}

BENCHMARKED_FUNCTION(principia_sin_cos) {
  auto const values = SinCos(x * Radian);
  // The nanobenchmark library wants the result to be a double, so we'll pay the
  // price of an extra `and` (1 cycle).
  return _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_sd(values.sin), _mm_set_sd(values.cos)));
}

}  // namespace _examples
}  // namespace nanobenchmarks
}  // namespace principia
