#include <emmintrin.h>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "nanobenchmarks/function_registry.hpp"  // 🧙 For BENCHMARK_FUNCTION etc.
#include "numerics/cbrt.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _examples {

using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_cbrt;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

BENCHMARKED_FUNCTION(twice) {
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

auto const c0 = Displacement<World>({0 * Metre, 0 * Metre, 1 * Metre});
auto const c1 = Velocity<World>(
    {0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second});
auto const c2 = Vector<Acceleration, World>({1 * Metre / Second / Second,
                                             0 * Metre / Second / Second,
                                             0 * Metre / Second / Second});
auto const c3 = Vector<Jerk, World>({1 * Metre / Second / Second / Second,
                                     0 * Metre / Second / Second / Second,
                                     0 * Metre / Second / Second / Second});
auto const c4 =
    Vector<Snap, World>({1 * Metre / Second / Second / Second / Second,
                         0 * Metre / Second / Second / Second / Second,
                         0 * Metre / Second / Second / Second / Second});
auto const c5 = Vector<Variation<Snap>, World>(
    {1 * Metre / Second / Second / Second / Second / Second,
     0 * Metre / Second / Second / Second / Second / Second,
     0 * Metre / Second / Second / Second / Second / Second});

BENCHMARKED_FUNCTION(polynomial1) {
  using P1A = PolynomialInMonomialBasis<Displacement<World>, Instant, 1>;

  P1A::Coefficients const coefficients({c0, c1});
  Instant const t0 = Instant() + 0.3 * Second;
  P1A const p(coefficients, t0, with_evaluator<Estrin>);
  auto const result = p(t0 + x * Second).coordinates();
  return _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_pd(result.x / Metre, result.y / Metre),
                 _mm_set_sd(result.z / Metre)));
}

BENCHMARKED_FUNCTION(polynomial2) {
  using P2A = PolynomialInMonomialBasis<Displacement<World>, Instant, 2>;

  P2A::Coefficients const coefficients({c0, c1, c2});
  Instant const t0 = Instant() + 0.3 * Second;
  P2A const p(coefficients, t0, with_evaluator<Estrin>);
  auto const result = p(t0 + x * Second).coordinates();
  return _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_pd(result.x / Metre, result.y / Metre),
                 _mm_set_sd(result.z / Metre)));
}

BENCHMARKED_FUNCTION(polynomial3) {
  using P3A = PolynomialInMonomialBasis<Displacement<World>, Instant, 3>;

  P3A::Coefficients const coefficients({c0, c1, c2, c3});
  Instant const t0 = Instant() + 0.3 * Second;
  P3A const p(coefficients, t0, with_evaluator<Estrin>);
  auto const result = p(t0 + x * Second).coordinates();
  return _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_pd(result.x / Metre, result.y / Metre),
                 _mm_set_sd(result.z / Metre)));
}

BENCHMARKED_FUNCTION(polynomial4) {
  using P4A = PolynomialInMonomialBasis<Displacement<World>, Instant, 4>;

  P4A::Coefficients const coefficients({c0, c1, c2, c3, c4});
  Instant const t0 = Instant() + 0.3 * Second;
  P4A const p(coefficients, t0, with_evaluator<Estrin>);
  auto const result = p(t0 + x * Second).coordinates();
  return _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_pd(result.x / Metre, result.y / Metre),
                 _mm_set_sd(result.z / Metre)));
}

BENCHMARKED_FUNCTION(polynomial5) {
  using P5A = PolynomialInMonomialBasis<Displacement<World>, Instant, 5>;

  P5A::Coefficients const coefficients({c0, c1, c2, c3, c4, c5});
  Instant const t0 = Instant() + 0.3 * Second;
  P5A const p(coefficients, t0, with_evaluator<Estrin>);
  auto const position = p(t0 + x * Second).coordinates();
  auto const velocity = p.EvaluateDerivative(t0 + x * Second).coordinates();
  return _mm_cvtsd_f64(_mm_and_pd(
      _mm_and_pd(_mm_set_pd(position.x / Metre, position.y / Metre),
                 _mm_set_pd(position.z / Metre, velocity.x / (Metre / Second))),
      _mm_set_pd(velocity.y / (Metre / Second),
                 velocity.z / (Metre / Second))));
}

}  // namespace _examples
}  // namespace nanobenchmarks
}  // namespace principia
