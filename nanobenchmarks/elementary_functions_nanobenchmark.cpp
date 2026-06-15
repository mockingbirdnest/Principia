#include <immintrin.h>

#include <cmath>

#include "core-math/cos.h"
#include "core-math/sin.h"
#include "nanobenchmarks/nanobenchmark.hpp"  // 🧙 For NANOBENCHMARK_*.
#include "numerics/angle_reduction.hpp"
#include "numerics/cbrt.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace nanobenchmarks {
namespace {

using namespace principia::numerics::_angle_reduction;
using namespace principia::numerics::_cbrt;
using namespace principia::numerics::_double_precision;
using namespace principia::nanobenchmarks::_nanobenchmark;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

NANOBENCHMARK_FUNCTION(Cbrt);

using namespace principia::numerics::_cbrt::internal;

NANOBENCHMARK_FUNCTION(method_3²ᴄZ5¹::Cbrt<Rounding::Faithful>);
NANOBENCHMARK_FUNCTION(method_3²ᴄZ5¹::Cbrt<Rounding::Correct>);
#if PRINCIPIA_COMPILER_MSVC
NANOBENCHMARK_FUNCTION(method_5²Z4¹FMA::Cbrt<Rounding::Faithful>);
NANOBENCHMARK_FUNCTION(method_5²Z4¹FMA::Cbrt<Rounding::Correct>);
#endif

NANOBENCHMARK(std_sin) {
  return std::sin(x);
}

NANOBENCHMARK(cr_sin) {
  return cr_sin(x);
}

NANOBENCHMARK(principia_sin) {
  return Sin(x * Radian);
}

NANOBENCHMARK(std_cos) {
  return std::cos(x);
}

NANOBENCHMARK(cr_cos) {
  return cr_cos(x);
}

NANOBENCHMARK(principia_cos) {
  return Cos(x * Radian);
}

NANOBENCHMARK(principia_sin_cos) {
  auto const values = SinCos(x * Radian);
  // The nanobenchmark library wants the result to be a double, so we'll pay the
  // price of an extra `and` (1 cycle).
  return _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_sd(values.sin), _mm_set_sd(values.cos)));
}

NANOBENCHMARK(payne_hanek) {
  DoublePrecision<Angle> x_reduced;
  std::int64_t quadrant;
  PayneHanekReduction<61>(x * Radian, x_reduced, quadrant);
  return x_reduced.value / Radian;
}

}  // namespace
}  // namespace nanobenchmarks
}  // namespace principia
