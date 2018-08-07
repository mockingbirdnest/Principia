
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

// Convenient overloads of the intrinsic functions.
// TODO(phl): Move to a central place?

__m128d _mm_mul_sd(double const left, double const right) {
  return _mm_mul_sd(_mm_set_sd(left), _mm_set_sd(right));
}

double _mm_mul_sd(double const left, __m128d const right) {
  return _mm_cvtsd_f64(_mm_mul_sd(_mm_set_sd(left), right));
}

__m128d _mm_sub_sd(__m128d const left, std::int64_t const right) {
  return _mm_sub_sd(left, _mm_cvtsi64_sd(left, right));
}

}  // namespace

void FastSinCos2π(double cycles, double& sin, double& cos) {
  // Argument reduction.  cycles_reduced is in [-1/8, 1/8] and quadrant goes
  // from 0 to 3, with 0 indicating the principal quadrant and the others
  // numbered in the trigonometric direction.  These quantities are computed by
  // extracting the integer and fractional parts of 4 * cycles.
  double cycles_reduced;
  std::int64_t four_cycles_integer;
#if PRINCIPIA_USE_SSE3_INTRINSICS
  __m128d const four_cycles = _mm_mul_sd(4.0, cycles);
  four_cycles_integer = _mm_cvtsd_si64(four_cycles);
  __m128d const four_cycles_fractional = _mm_sub_sd(four_cycles,
                                                    four_cycles_integer);
  cycles_reduced = _mm_mul_sd(0.25, four_cycles_fractional);
#else
  double const four_cycles = 4.0 * cycles;
  four_cycles_integer = std::nearbyint(four_cycles);
  cycles_reduced = 0.25 * (four_cycles - four_cycles_integer);
#endif
  std::int64_t const quadrant = four_cycles_integer & 0b11;

  double const cycles_reduced² = cycles_reduced * cycles_reduced;
  double const s = sin_polynomial.Evaluate(cycles_reduced²) * cycles_reduced;
  double const c = 1.0 +
                   cos_polynomial.Evaluate(cycles_reduced²) * cycles_reduced²;

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
