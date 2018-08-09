
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
                                     /*Evaluator=*/EstrinEvaluator>;

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

struct Decomposition {
  std::int64_t integer_part;
  double fractional_part;
};

Decomposition Decompose(double const x) {
  Decomposition decomposition;
#if PRINCIPIA_USE_SSE3_INTRINSICS
  __m128d const x_128d = _mm_set_sd(x);
  decomposition.integer_part = _mm_cvtsd_si64(x_128d);
  decomposition.fractional_part = _mm_cvtsd_f64(
      _mm_sub_sd(x_128d,
                 _mm_cvtsi64_sd(__m128d{}, decomposition.integer_part)));
#else
  decomposition.integer_part = std::nearbyint(x);
  decomposition.fractional_part = x - decomposition.integer_part;
#endif
  return decomposition;
}

}  // namespace

void FastSinCos2π(double cycles, double& sin, double& cos) {
  // Argument reduction.  cycles_reduced is in [-1/8, 1/8] and quadrant goes
  // from 0 to 3, with 0 indicating the principal quadrant and the others
  // numbered in the trigonometric direction.  These quantities are computed by
  // extracting the integer and fractional parts of 4 * cycles.
  Decomposition const decomposition = Decompose(4.0 * cycles);
  double const cycles_reduced = 0.25 * decomposition.fractional_part;
  std::int64_t const quadrant = decomposition.integer_part & 0b11;

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
