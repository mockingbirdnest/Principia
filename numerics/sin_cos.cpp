#pragma once

#include "numerics/sin_cos.hpp"

#include <pmmintrin.h>

#include "numerics/accurate_tables.mathematica.h"
#include "numerics/double_precision.hpp"
#include "numerics/fma.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

using namespace principia::numerics::_accurate_tables;
using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_fma;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::quantities::_elementary_functions;

using Argument = double;
using Value = double;

constexpr Argument sin_near_zero_cutoff = 1.0 / 1024.0;
constexpr Argument table_spacing_reciprocal = 512.0;

// π / 2 split so that the high half has 18 zeros at the end of its mantissa.
constexpr std::int64_t π_over_2_zeroes = 18;
constexpr Argument π_over_2_threshold = π / 2 * (1 << π_over_2_zeroes);
constexpr Argument π_over_2_high = 0x1.921fb'5444p0;
constexpr Argument π_over_2_low = 0x2.d1846'9898c'c5170'1b839p-40;

template<FMAPolicy fma_policy>
using Polynomial1 = HornerEvaluator<Value, Argument, 1, fma_policy>;

namespace masks {
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
}  // namespace masks

template<FMAPolicy fma_policy>
double FusedMultiplyAdd(double const a, double const b, double const c) {
  if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
      (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using quantities::_elementary_functions::FusedMultiplyAdd;
    return FusedMultiplyAdd(a, b, c);
  } else {
    return a * b + c;
  }
}

void Reduce(Argument const x,
            DoublePrecision<Argument>& x_reduced,
            std::int64_t& quadrant) {
  if (x < π / 4 && x > -π / 4) {
    x_reduced.value = x;
    x_reduced.error = 0;
    quadrant = 0;
  } else if (x <= π_over_2_threshold && x >= -π_over_2_threshold) {
    std::int64_t const n = _mm_cvtsd_si64(_mm_set_sd(x * (2 / π)));
    double const n_double = static_cast<double>(n);
    Argument const value = x - n_double * π_over_2_high;
    Argument const error = n_double * π_over_2_low;
    x_reduced = QuickTwoDifference(value, error);
    // TODO(phl): Check for value too small.
    quadrant = n & 0b11;
  }
  // TODO(phl): Fallback to a more precise algorithm.
}

// TODO(phl): Take the perturbation into account in the polynomials.

template<FMAPolicy fma_policy>
Value SinPolynomial(Argument const x) {
  // Absolute error better than 84.8 bits over an interval of radius 1/1024.
  return Polynomial1<fma_policy>::Evaluate(
      {-0x1.5555555555555'555p-3, 0x1.111110B24ACB5'617p-7}, x);
}

template<FMAPolicy fma_policy>
Value SinPolynomialNearZero(Argument const x) {
  // Relative error better than 74.5 bits over an interval of radius 1/1024.
  return Polynomial1<fma_policy>::Evaluate(
      {-0x1.5555555555555'555p-3, 0x1.111110B40E889'1076p-7}, x);
}

template<FMAPolicy fma_policy>
Value CosPolynomial(Argument const x) {
  // Absolute error better than 72.4 bits over an interval of radius 1/1024.
  return Polynomial1<fma_policy>::Evaluate(
      {-0x1.FFFFFFFFFFFFF'000p-2, 0x1.555554B290E69'14113p-5}, x);
}

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value SinImplementation(DoublePrecision<Argument> const argument) {
  auto const& x = argument.value;
  auto const& e = argument.error;
  double const abs_x = std::abs(x);
  if (abs_x < sin_near_zero_cutoff) {
    double const x² = x * x;
    double const x³ = x² * x;
    return x + FusedMultiplyAdd<fma_policy>(
                   x³, SinPolynomialNearZero<fma_policy>(x²), e);
  } else {
    __m128d const sign = _mm_and_pd(masks::sign_bit, _mm_set_sd(x));
    auto const i = _mm_cvtsd_si64(_mm_set_sd(abs_x * table_spacing_reciprocal));
    auto const& accurate_values = SinCosAccurateTable[i];
    double const x₀ =
        _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.x), sign));
    double const sin_x₀ =
        _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.sin_x), sign));
    double const& cos_x₀ = accurate_values.cos_x;
  // [GB91] incorporates `e` in the computation of `h`.  However, `e` can be
  // very small and `h` can be as large as 1/512.  We note that the only
  // influence of `e` is on the computation of `cos_x₀ * e`, so we incorporate
  // `e` with the `SinPolynomial` below.
    double const h = x - x₀;

    DoublePrecision<double> const sin_x₀_plus_h_cos_x₀ =
        TwoProductAdd<fma_policy>(cos_x₀, h, sin_x₀);
    double const h² = h * h;
    double const h³ = h² * h;
    return sin_x₀_plus_h_cos_x₀.value +
           ((sin_x₀ * h² * CosPolynomial<fma_policy>(h²) +
             cos_x₀ * FusedMultiplyAdd<fma_policy>(
                          h³, SinPolynomial<fma_policy>(h²), e)) +
            sin_x₀_plus_h_cos_x₀.error);
  }
}

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value CosImplementation(DoublePrecision<Argument> const argument) {
  auto const& x = argument.value;
  auto const& e = argument.error;
  double const abs_x = std::abs(x);
  __m128d const sign = _mm_and_pd(masks::sign_bit, _mm_set_sd(x));
  auto const i = _mm_cvtsd_si64(_mm_set_sd(abs_x * table_spacing_reciprocal));
  auto const& accurate_values = SinCosAccurateTable[i];
  double const x₀ =
      _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.x), sign));
  double const sin_x₀ =
      _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.sin_x), sign));
  double const& cos_x₀ = accurate_values.cos_x;
  // [GB91] incorporates `e` in the computation of `h`.  However, `e` can be
  // very small and `h` can be as large as 1/512.  We note that the only
  // influence of `e` is on the computation of `sin_x₀ * e`, so we incorporate
  // `e` with the `SinPolynomial` below.
  double const h = x - x₀;

  DoublePrecision<double> const cos_x₀_minus_h_sin_x₀ =
      TwoProductNegatedAdd<fma_policy>(sin_x₀, h, cos_x₀);
  double const h² = h * h;
  double const h³ = h² * h;
  return cos_x₀_minus_h_sin_x₀.value +
         ((cos_x₀ * h² * CosPolynomial<fma_policy>(h²) -
           sin_x₀ * FusedMultiplyAdd<fma_policy>(
                        h³, SinPolynomial<fma_policy>(h²), e)) +
          cos_x₀_minus_h_sin_x₀.error);
}

#if PRINCIPIA_INLINE_SIN_COS
inline
#endif
Value __cdecl Sin(Argument const x) {
  DoublePrecision<Argument> x_reduced;
  std::int64_t quadrant;
  Reduce(x, x_reduced, quadrant);
  double value;
  if (UseHardwareFMA) {
    if (quadrant & 0b1) {
      value = CosImplementation<FMAPolicy::Force>(x_reduced);
    } else {
      value = SinImplementation<FMAPolicy::Force>(x_reduced);
    }
  } else {
    if (quadrant & 0b1) {
      value = CosImplementation<FMAPolicy::Disallow>(x_reduced);
    } else {
      value = SinImplementation<FMAPolicy::Disallow>(x_reduced);
    }
  }
  if (quadrant & 0b10) {
    return -value;
  } else {
    return value;
  }
}

#if PRINCIPIA_INLINE_SIN_COS
inline
#endif
Value __cdecl Cos(Argument const x) {
  DoublePrecision<Argument> x_reduced;
  std::int64_t quadrant;
  Reduce(x, x_reduced, quadrant);
  double value;
  if (UseHardwareFMA) {
    if (quadrant & 0b1) {
      value = SinImplementation<FMAPolicy::Force>(x_reduced);
    } else {
      value = CosImplementation<FMAPolicy::Force>(x_reduced);
    }
  } else {
    if (quadrant & 0b1) {
      value = SinImplementation<FMAPolicy::Disallow>(x_reduced);
    } else {
      value = CosImplementation<FMAPolicy::Disallow>(x_reduced);
    }
  }
  if (quadrant == 1 || quadrant == 2) {
    return -value;
  } else {
    return value;
  }
}

}  // namespace internal
}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
