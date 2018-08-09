
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
// The coefficients are scaled to accept arguments scaled to [0, 1/2].

// The sine polynomial uses a custom evaluation, so the individual coefficients
// are named.  The polynomial for Sin(2 π x/4) is s₁ x + s₃ x³ + s₅ x⁵.
double const s₁ = 6.28315387593158874093559349802 / 4;
double const s₃ = -41.3255673715186216778612605095 / (4 * 16);
double const s₅ = 79.5314110676979262924240784281 / (4 * 16 * 16);

// The cosine polynomial for 16 x² ⟼ (Cos(2 π x) - 1)/(16 x²) is turned into a
// polynomial for 16 x² ⟼ Cos(2 π x) by multiplying by the argument and adding
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
  // Argument reduction.  cycles_reduced is in [-1/2, 1/2] and quadrant goes
  // from 0 to 3, with 0 indicating the principal quadrant and the others
  // numbered in the trigonometric direction.  These quantities are computed by
  // extracting the integer and fractional parts of 4 * cycles.
  Decomposition const decomposition = Decompose(4.0 * cycles);
  // TODO(egg): cycles_reduced is no longer an appropriate name.
  // This now counts a fraction of a right angle/quarter cycle/quadrant.
  double const cycles_reduced = decomposition.fractional_part;
  double const cycles_reduced² = cycles_reduced * cycles_reduced;
  double const cycles_reduced³ = cycles_reduced² * cycles_reduced;
  std::int64_t const quadrant = decomposition.integer_part & 0b11;

  // The custom evaluation, compared to Estrin followed by multiplication by the
  // argument, i.e., x * (s₁ + s₃ * x² + s₅ * (x² * x²)), avoids having a
  // multiplication by x at the end.
  double const s =
      s₁ * cycles_reduced + (s₃ + s₅ * cycles_reduced²) * cycles_reduced³;
  // TODO(egg): Considering the above, I'm not sure this is more readable than
  // just spelling out the Estrin evaluation.
  double const c = cos_polynomial.Evaluate(cycles_reduced²);

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
