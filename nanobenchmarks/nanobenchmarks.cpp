#include <emmintrin.h>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "nanobenchmarks/function_registry.hpp"  // 🧙 For NANOBENCHMARK_*.
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
using namespace principia::nanobenchmarks::_function_registry;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

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

NANOBENCHMARK(principia_sin) {
  return Sin(x * Radian);
}

NANOBENCHMARK(std_cos) {
  return std::cos(x);
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

class PolynomialNanobenchmark : public Fixture {
 protected:
  PolynomialNanobenchmark()
      : c0_({0 * Metre, 0 * Metre, 1 * Metre}),
        c1_({0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second}),
        c2_({1 * Metre / Second / Second,
             0 * Metre / Second / Second,
             0 * Metre / Second / Second}),
        c3_({1 * Metre / Second / Second / Second,
             0 * Metre / Second / Second / Second,
             0 * Metre / Second / Second / Second}),
        c4_({1 * Metre / Second / Second / Second / Second,
             0 * Metre / Second / Second / Second / Second,
             0 * Metre / Second / Second / Second / Second}),
        c5_({1 * Metre / Second / Second / Second / Second / Second,
             0 * Metre / Second / Second / Second / Second / Second,
             0 * Metre / Second / Second / Second / Second / Second}) {};

  Displacement<World> const c0_;
  Velocity<World> const c1_;
  Vector<Acceleration, World> const c2_;
  Vector<Jerk, World> const c3_;
  Vector<Snap, World> const c4_;
  Vector<Variation<Snap>, World> const c5_;
};

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree1) {
  using P1A = PolynomialInMonomialBasis<Displacement<World>, Instant, 1>;

  P1A::Coefficients const coefficients({c0_, c1_});
  Instant const t0 = Instant() + 0.3 * Second;
  P1A const p(coefficients, t0, with_evaluator<Estrin>);
  return p(t0 + x * Second).Norm²() / Metre / Metre;
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree2) {
  using P2A = PolynomialInMonomialBasis<Displacement<World>, Instant, 2>;

  P2A::Coefficients const coefficients({c0_, c1_, c2_});
  Instant const t0 = Instant() + 0.3 * Second;
  P2A const p(coefficients, t0, with_evaluator<Estrin>);
  return p(t0 + x * Second).Norm²() / Metre / Metre;
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree3) {
  using P3A = PolynomialInMonomialBasis<Displacement<World>, Instant, 3>;

  P3A::Coefficients const coefficients({c0_, c1_, c2_, c3_});
  Instant const t0 = Instant() + 0.3 * Second;
  P3A const p(coefficients, t0, with_evaluator<Estrin>);
  return p(t0 + x * Second).Norm²() / Metre / Metre;
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree4) {
  using P4A = PolynomialInMonomialBasis<Displacement<World>, Instant, 4>;

  P4A::Coefficients const coefficients({c0_, c1_, c2_, c3_, c4_});
  Instant const t0 = Instant() + 0.3 * Second;
  P4A const p(coefficients, t0, with_evaluator<Estrin>);
  return p(t0 + x * Second).Norm²() / Metre / Metre;
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree5) {
  using P5A = PolynomialInMonomialBasis<Displacement<World>, Instant, 5>;

  P5A::Coefficients const coefficients({c0_, c1_, c2_, c3_, c4_, c5_});
  Instant const t0 = Instant() + 0.3 * Second;
  P5A const p(coefficients, t0, with_evaluator<Estrin>);
  return p(t0 + x * Second).Norm²() / Metre / Metre;
}

}  // namespace _examples
}  // namespace nanobenchmarks
}  // namespace principia
