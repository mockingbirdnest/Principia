
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

// 2nd-degree polynomials that minimize the absolute error on sin and cos over
// the interval [0, 1/8].  The minimization algorithm is run on
// Sin(2 π √x)/√x and (Cos(2 π √x) - 1)/x to ensure that the functions have
// the right behavior near 0 and the proper parity.  Because of extra
// oscillations, the lower bounds of the minimization intervals are 1/36 and
// 1/24 respectively.  This is where the maximum error is found.
// The coefficients are scaled to accept arguments expressed in right angles.

// The sine polynomial uses a custom evaluation, so the individual coefficients
// are named.  The polynomial for Sin(2 π y/4) is s₁ y + s₃ y³ + s₅ y⁵.
double const s₁ = 6.28315387593158874093559349802 / 4;
double const s₃ = -41.3255673715186216778612605095 / (4 * 16);
double const s₅ = 79.5314110676979262924240784281 / (4 * 16 * 16);

// The cosine polynomial for 16 z² ⟼ (Cos(2 π z) - 1)/(16 z²) is turned into a
// polynomial for 16 z² ⟼ Cos(2 π z) by multiplying by the argument and adding
// 1, i.e., prepending 1 to the list of coefficients.
P3 cos_polynomial(P3::Coefficients{
    1.0,
    -19.7391672615468690589481752820 / 16,
    64.9232282990046449731568966307 / (16 * 16),
    -83.6659064641344641438100039739 / (16 * 16 * 16)});

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

void FastSinCos2π(double const cycles, double& sin, double& cos) {
  // Argument reduction.
  // - quadrant goes from 0 to 3, with 0 indicating the principal quadrant and
  //   the others numbered in the trigonometric direction;
  // - y is in [-1/2, 1/2], corresponding to reduced angles of [-π/4, π/4] in
  //   the quadrant.
  // These quantities are computed by extracting the integer and fractional
  // parts of 4 * cycles.
  Decomposition const decomposition = Decompose(4.0 * cycles);
  double const y = decomposition.fractional_part;
  double const y² = y * y;
  double const y³ = y² * y;
  std::int64_t const quadrant = decomposition.integer_part & 0b11;

  // The custom evaluation, compared to Estrin followed by multiplication by the
  // argument, i.e., y * (s₁ + s₃ * y² + s₅ * (y² * y²)), avoids having a
  // multiplication by y at the end.
  double const s = s₁ * y + (s₃ + s₅ * y²) * y³;
  double const c = cos_polynomial.Evaluate(y²);

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
