
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
  std::int64_t quadrant;  // 0..3
#if PRINCIPIA_USE_SSE3_INTRINSICS
  __m128d const cycles_128d = _mm_load1_pd(&cycles);
  __m128d const two_cycles_128d = _mm_add_sd(cycles_128d, cycles_128d);
  __m128d const four_cycles_128d = _mm_add_sd(two_cycles_128d, two_cycles_128d);
  __int64 const four_cycles_64 = _mm_cvtsd_si64(four_cycles_128d);
  quadrant = four_cycles_64 & 0b11;
  __m128d const four_cycles_integer_128d =
      _mm_cvtsi64_sd(four_cycles_128d, four_cycles_64);
  __m128d const four_cycles_fractional_128d =
      _mm_sub_sd(four_cycles_128d, four_cycles_integer_128d);
  __m128d const one_fourth = _mm_set1_pd(0.25);
  __m128d const cycles_fractional_128d =
      _mm_mul_sd(one_fourth, four_cycles_fractional_128d);
  cycles_fractional = _mm_cvtsd_f64(cycles_fractional_128d);
#else
  double const cycles_integer = std::nearbyint(cycles);
  cycles_fractional = cycles - cycles_integer;
#endif

  double const cycles_fractional² = cycles_fractional * cycles_fractional;
  double const s =
      sin_polynomial.Evaluate(cycles_fractional²) * cycles_fractional;
  double const c =
      1.0 + cos_polynomial.Evaluate(cycles_fractional²) * cycles_fractional²;

  //LOG(ERROR)<<quadrant<<" "<<cycles_fractional;
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
