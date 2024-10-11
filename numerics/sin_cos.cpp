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
static const __m128d exponent_bits =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x7ff0'0000'0000'0000));
static const __m128d mantissa_bits =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x000f'ffff'ffff'ffff));
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

int alerts = 0;

double DoTheThing(double const x, double const dx) {
  double const value = x + dx;
  double const error = (value - x) - dx;
  __m128i const value_exponent_128i =
      _mm_castpd_si128(_mm_and_pd(masks::exponent_bits, _mm_set_sd(value)));
  double const value_exponent =
      _mm_cvtsd_f64(_mm_castsi128_pd(value_exponent_128i));
  __m128i const error_128i =
      _mm_castpd_si128(_mm_andnot_pd(masks::sign_bit, _mm_set_sd(error)));
  double const normalized_error = _mm_cvtsd_f64(
      _mm_castsi128_pd(_mm_sub_epi64(error_128i, value_exponent_128i)));
  if (normalized_error < -0x1.ffffp971 && normalized_error > -0x1.00008p972) {
    //LOG(ERROR) << std::setprecision(25) << x << " " << std::hexfloat
    //           << value<<" "<<error<< " "<< normalized_error;
    ++alerts;
  }
  return value;
}

void Reduce(Argument const x,
            DoublePrecision<Argument>& x_reduced,
            std::int64_t& quadrant) {
  if (x < π / 4 && x > -π / 4) {
    x_reduced.value = x;
    x_reduced.error = 0;
    quadrant = 0;
  } else if (x <= π_over_2_threshold && x >= -π_over_2_threshold) {
    // We are not very sensitive to rounding errors in this expression, because
    // in the worst case it could cause the reduced angle to jump from the
    // vicinity of π / 4 to the vicinity of -π / 4 with appropriate adjustment
    // of the quadrant.
    __m128d const n_128d = _mm_round_sd(
        _mm_setzero_pd(), _mm_set_sd(x * (2 / π)), _MM_FROUND_RINT);
    double const n_double = _mm_cvtsd_f64(n_128d);
    std::int64_t const n = _mm_cvtsd_si64(n_128d);
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
    double const x³_term = FusedMultiplyAdd<fma_policy>(
                   x³, SinPolynomialNearZero<fma_policy>(x²), e);
    return DoTheThing(x, x³_term);
  } else {
    __m128d const sign = _mm_and_pd(masks::sign_bit, _mm_set_sd(x));
    auto const i = _mm_cvtsd_si64(_mm_set_sd(abs_x * table_spacing_reciprocal));
    auto const& accurate_values = SinCosAccurateTable[i];
    double const x₀ =
        _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.x), sign));
    double const sin_x₀ =
        _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.sin_x), sign));
    double const& cos_x₀ = accurate_values.cos_x;
    // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
    // don't overlap and in the first interval `x` and `h` may be of the same
    // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
    // later in the computation.
    double const h = x - x₀;

    DoublePrecision<double> const sin_x₀_plus_h_cos_x₀ =
        TwoProductAdd<fma_policy>(cos_x₀, h, sin_x₀);
    double const h² = h * (h + 2 * e);
    double const h³ = h² * h;
    double const polynomial_term =
        ((sin_x₀ * h² * CosPolynomial<fma_policy>(h²) +
          cos_x₀ * FusedMultiplyAdd<fma_policy>(
                       h³, SinPolynomial<fma_policy>(h²), e)) +
         sin_x₀_plus_h_cos_x₀.error);
    return DoTheThing(sin_x₀_plus_h_cos_x₀.value, polynomial_term);
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
  // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
  // don't overlap and in the first interval `x` and `h` may be of the same
  // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
  // later in the computation.
  double const h = x - x₀;

  DoublePrecision<double> const cos_x₀_minus_h_sin_x₀ =
      TwoProductNegatedAdd<fma_policy>(sin_x₀, h, cos_x₀);
  double const h² = h * (h + 2 * e);
  double const h³ = h² * h;
  double const polynomial_term =
      ((cos_x₀ * h² * CosPolynomial<fma_policy>(h²) -
        sin_x₀ * FusedMultiplyAdd<fma_policy>(
                     h³, SinPolynomial<fma_policy>(h²), e)) +
       cos_x₀_minus_h_sin_x₀.error);
  return DoTheThing(cos_x₀_minus_h_sin_x₀.value, polynomial_term);
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
