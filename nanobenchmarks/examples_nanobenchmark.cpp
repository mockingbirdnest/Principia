#include <emmintrin.h>

#include "nanobenchmarks/function_registry.hpp"  // 🧙 For NANOBENCHMARK_*.

namespace principia {
namespace nanobenchmarks {
namespace _examples {

using namespace principia::nanobenchmarks::_function_registry;

NANOBENCHMARK(twice) {
  return 2 * x;
}

NANOBENCHMARK(thrice) {
  return 3 * x;
}

NANOBENCHMARK(inc) {
  return x + 1;
}

NANOBENCHMARK(multiply_4_times) {
  return x * x * x * x * x;
}

NANOBENCHMARK(add_16_times) {
  return x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x;
}

NANOBENCHMARK(square_root) {
  __m128d x_0 = _mm_set_sd(x);
  return _mm_cvtsd_f64(_mm_sqrt_sd(x_0, x_0));
}

NANOBENCHMARK(sqrt_sqrt) {
  __m128d x_0 = _mm_set_sd(x);
  x_0 = _mm_sqrt_sd(x_0, x_0);
  return _mm_cvtsd_f64(_mm_sqrt_sd(x_0, x_0));
}

NANOBENCHMARK(square_root_division) {
  __m128d x_0 = _mm_set_sd(x);
  return _mm_cvtsd_f64(_mm_div_sd(x_0, _mm_sqrt_sd(x_0, x_0)));
}

}  // namespace _examples
}  // namespace nanobenchmarks
}  // namespace principia
