
#pragma once

#include "numerics/fast_sin_cos_2π.hpp"

#include <pmmintrin.h>

#include "base/macros.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"

namespace principia {
namespace numerics {

namespace {

using P3 = PolynomialInMonomialBasis</*Value=*/double,
                                     /*Argument=*/double,
                                     /*degree=*/3,
                                     /*Evaluator=*/EstrinEvaluator>;

// 3rd-degree polynomials that minimize the absolute error on sin and cos over
// the interval [0, 1/4].  The minimization algorithm is run on
// Sin(2 π √x)/√x and (Cos(2 π √x) - 1)/x to ensure that the functions have
// the right behavior near 0 and the proper parity.  Because of extra
// oscillations, the lower bounds of the minimization intervals are 1/23 and
// 1/15 respectively.  This is where the maximum error is found.
P3 sin_polynomial(P3::Coefficients{6.28316404405113818577981340506,
                                   -41.3371423477858688509416864345,
                                   81.3407682603799599938651480917,
                                   -70.9934281315300026308830925659});
P3 cos_polynomial(P3::Coefficients{-19.7391820689166533085010275514,
                                   64.9352211775039525420259682190,
                                   -85.2540004035261113433497714557,
                                   56.3405940928237075549636782571});

}  // namespace

void FastSinCos2π(double cycles, double& sin, double& cos) {
  // Argument reduction.  Since the argument is in cycles, we just drop the
  // integer part.
  double cycles_fractional;
#if PRINCIPIA_USE_SSE3_INTRINSICS
  __m128d const cycles_128d = _mm_load1_pd(&cycles);
  __int64 const cycles_64 = _mm_cvtsd_si64(cycles_128d);
  __m128d const cycles_integer_128d = _mm_cvtsi64_sd(cycles_128d, cycles_64);
  __m128d const cycles_fractional_128d = _mm_sub_sd(cycles_128d,
                                                    cycles_integer_128d);
  cycles_fractional = _mm_cvtsd_f64(cycles_fractional_128d);
#else
  double const cycles_integer = std::nearbyint(cycles);
  cycles_fractional = cycles - cycles_integer;
#endif

  double sign = 1.0;
  if (cycles_fractional > 0.25) {
    cycles_fractional -= 0.5;
    sign = -1.0;
  } else if (cycles_fractional < -0.25) {
    cycles_fractional += 0.5;
    sign = -1.0;
  }

  double const cycles_fractional² = cycles_fractional * cycles_fractional;
  sin = sin_polynomial.Evaluate(cycles_fractional²) *
        (sign * cycles_fractional);
  cos = sign + cos_polynomial.Evaluate(cycles_fractional²) *
               (sign * cycles_fractional²);
}

}  // namespace numerics
}  // namespace principia
