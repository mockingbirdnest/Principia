
#pragma once

#include "numerics/fast_sin_cos_cycle.hpp"

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
P3 sin_polynomial(P3::Coefficients{6.28316404405113818577981340506,
                                   -41.3371423477858688509416864345,
                                   81.3407682603799599938651480917,
                                   -70.9934281315300026308830925659});
P3 cos_polynomial(P3::Coefficients{-19.7391820689166533085010275514,
                                   64.9352211775039525420259682190,
                                   -85.2540004035261113433497714557,
                                   56.3405940928237075549636782571});

}  // namespace

void FastSinCosCycle(double x, double& sin, double& cos) {
  // Argument reduction.  Since the unit of x is the cycle, we just drop the
  // integer part.
  double x_fractional;
#if PRINCIPIA_USE_SSE3_INTRINSICS
  __m128d const x_128d = _mm_load1_pd(&x);
  // Note that _mm_cvtpd has lower latency than _mm_cvtsd.
  __m128i const x_128i = _mm_cvtpd_epi32(x_128d);
  __m128d const x_integer_128d = _mm_cvtepi32_pd(x_128i);
  __m128d const x_fractional_128d = _mm_sub_pd(x_128d, x_integer_128d);
  x_fractional = _mm_cvtsd_f64(x_fractional_128d);
#else
#error "Argument reduction not implemented for non-SSE3 processors"
#endif

  double sign = 1.0;
  if (x_fractional > 0.25) {
    x_fractional -= 0.5;
    sign = -1.0;
  } else if (x_fractional < -0.25) {
    x_fractional += 0.5;
    sign = -1.0;
  }

  double const x_fractional² = x_fractional * x_fractional;
  sin = sin_polynomial.Evaluate(x_fractional²) * (sign * x_fractional);
  cos = sign + cos_polynomial.Evaluate(x_fractional²) * (sign * x_fractional²);
}

}  // namespace numerics
}  // namespace principia
