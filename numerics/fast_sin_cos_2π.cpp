
#pragma once

#include "numerics/fast_sin_cos_2π.hpp"

#include <pmmintrin.h>

#include "base/macros.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"

namespace principia {
namespace numerics {

namespace {

using P2 = PolynomialInMonomialBasis</*Value=*/double,
                                     /*Argument=*/double,
                                     /*degree=*/2,
                                     /*Evaluator=*/HornerEvaluator>;

// 2nd-degree polynomials that minimize the absolute error on sin and cos over
// the interval [0, 1/8].  The minimization algorithm is run on
// Sin(2 π √x)/√x and (Cos(2 π √x) - 1)/x to ensure that the functions have
// the right behavior near 0 and the proper parity.  Because of extra
// oscillations, the lower bounds of the minimization intervals are 1/36 and
// 1/24 respectively.  This is where the maximum error is found.
P2 sin_polynomial(P2::Coefficients{6.28315387593158874093559349802,
                                   -41.3255673715186216778612605095,
                                   79.5314110676979262924240784281});
P2 cos_polynomial(P2::Coefficients{-19.7391672615468690589481752820,
                                   64.9232282990046449731568966307,
                                   -83.6659064641344641438100039739});

}  // namespace

void FastSinCos2π(double cycles, double& sin, double& cos) {
  // Argument reduction.  Since the argument is in cycles, we just drop the
  // integer part.
  //TODO(phl):comment.
  double cycles_fractional;
  std::int64_t quadrant;  // 0..3
#if PRINCIPIA_USE_SSE3_INTRINSICS
  __m128d const cycles_128d = _mm_set_sd(cycles);
  __m128d const four = _mm_set_sd(4.0);
  __m128d const four_cycles_128d = _mm_mul_sd(four, cycles_128d);
  __int64 const four_cycles_64 = _mm_cvtsd_si64(four_cycles_128d);
  quadrant = four_cycles_64 & 0b11;
  __m128d const four_cycles_integer_128d =
      _mm_cvtsi64_sd(four_cycles_128d, four_cycles_64);
  __m128d const four_cycles_fractional_128d =
      _mm_sub_sd(four_cycles_128d, four_cycles_integer_128d);
  __m128d const one_fourth = _mm_set_sd(0.25);
  __m128d const cycles_fractional_128d =
      _mm_mul_sd(one_fourth, four_cycles_fractional_128d);
  cycles_fractional = _mm_cvtsd_f64(cycles_fractional_128d);
#else
  double const cycles_integer = std::nearbyint(cycles);
  cycles_fractional = cycles - cycles_integer;
#endif
  //TODO(phl):debug.

  double const cycles_fractional² = cycles_fractional * cycles_fractional;
  double const s =
      sin_polynomial.Evaluate(cycles_fractional²) * cycles_fractional;
  double const c =
      1.0 + cos_polynomial.Evaluate(cycles_fractional²) * cycles_fractional²;

  switch (quadrant) {
    case 0:
      sin = s;
      cos = c;
      break;
    case 1:
      sin = c;
      cos = -s;
      break;
    case 2:
      sin = -s;
      cos = -c;
      break;
    case 3:
      sin = -c;
      cos = s;
      break;
  }
}

}  // namespace numerics
}  // namespace principia
