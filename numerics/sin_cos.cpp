#pragma once

#include "numerics/sin_cos.hpp"

#include <pmmintrin.h>

#include <limits>

#include "core-math/cos.h"
#include "core-math/sin.h"
#include "numerics/accurate_tables.mathematica.h"
#include "numerics/double_precision.hpp"
#include "numerics/fma.hpp"
#include "numerics/osaca.hpp"  // 🧙 For OSACA_*.
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

#define OSACA_ANALYSED_FUNCTION Cos
#define OSACA_ANALYSED_FUNCTION_NAMESPACE
#define OSACA_ANALYSED_FUNCTION_TEMPLATE_PARAMETERS
#define UNDER_OSACA_HYPOTHESES(expression)                                   \
  [&] {                                                                      \
    constexpr bool UseHardwareFMA = true;                                    \
    constexpr double θ = 3;                                                  \
    /* From argument reduction. */                                           \
    constexpr double abs_θ = θ > 0 ? θ : -θ;                                 \
    constexpr std::int64_t n = static_cast<std::int64_t>(θ * (2 / π) + 0.5); \
    constexpr double reduction_value = θ - n * π_over_2_high;                \
    constexpr double reduction_error = n * π_over_2_low;                     \
    /* Used to determine whether a better argument reduction is needed. */   \
    constexpr DoublePrecision<double> θ_reduced =                            \
        QuickTwoDifference(reduction_value, reduction_error);                \
    /* Used in Sin to detect the near-0 case. */                             \
    constexpr double abs_x =                                                 \
        θ_reduced.value > 0 ? θ_reduced.value : -θ_reduced.value;            \
    /* Used throughout the top-level functions. */                           \
    constexpr std::int64_t quadrant = n & 0b11;                              \
    /* Used in DetectDangerousRounding. */                                   \
    constexpr double normalized_error = 0;                                   \
    /* Not NaN is the only part that matters; used at the end of the    */   \
    /* top-level functions to determine whether to call the slow path.  */   \
    constexpr double value = 1;                                              \
    return expression;                                                       \
  }()

using Argument = double;
using Value = double;

constexpr std::int64_t table_spacing_bits = 9;
constexpr Argument table_spacing_reciprocal = 1 << table_spacing_bits;
constexpr Argument table_spacing = 1.0 / table_spacing_reciprocal;
constexpr Argument sin_near_zero_cutoff = table_spacing / 2.0;

// π / 2 split so that the high half has 18 zeros at the end of its mantissa.
constexpr std::int64_t π_over_2_zeroes = 18;
constexpr Argument π_over_2_threshold = π / 2 * (1 << π_over_2_zeroes);
constexpr Argument π_over_2_high = 0x1.921fb'5444p0;
constexpr Argument π_over_2_low = 0x2.d1846'9898c'c5170'1b839p-40;
constexpr Argument mantissa_reduce_shifter =
    1LL << (std::numeric_limits<double>::digits - 1);

template<FMAPolicy fma_policy>
using Polynomial1 = HornerEvaluator<Value, Argument, 1, fma_policy>;

namespace masks {
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d exponent_bits =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x7ff0'0000'0000'0000));
static const __m128d mantissa_bits =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x000f'ffff'ffff'ffff));
static const __m128d mantissa_index_bits =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x0000'0000'0000'01ff));
}  // namespace masks

inline std::int64_t AccurateTableIndex(double const abs_x) {
  // This function computes the index in the accurate table:
  // 1. A suitable (large) power of 2 is added to the argument so that the last
  //    bit of the mantissa of the result corresponds to units of 1/512 of the
  //    argument.  As part of this addition, an interval of radius 1/1024 around
  //    an integral multiple of 1/512 is correctly rounded to that integral
  //    multiple.
  // 2. An `and` operation is used to only retain the last 9 bits of the
  //    mantissa.
  // 3. The result is interpreted as an integer and returned as the index.
  return _mm_cvtsi128_si64(_mm_castpd_si128(_mm_and_pd(
      masks::mantissa_index_bits,
      _mm_set_sd(abs_x + (1LL << (std::numeric_limits<double>::digits -
                                  table_spacing_bits - 1))))));
}

template<FMAPolicy fma_policy>
double FusedMultiplyAdd(double const a, double const b, double const c) {
  static_assert(fma_policy != FMAPolicy::Auto);
  if constexpr (fma_policy == FMAPolicy::Force) {
    using quantities::_elementary_functions::FusedMultiplyAdd;
    return FusedMultiplyAdd(a, b, c);
  } else {
    return a * b + c;
  }
}

template<FMAPolicy fma_policy>
double FusedNegatedMultiplyAdd(double const a, double const b, double const c) {
  static_assert(fma_policy != FMAPolicy::Auto);
  if constexpr (fma_policy == FMAPolicy::Force) {
    using quantities::_elementary_functions::FusedNegatedMultiplyAdd;
    return FusedNegatedMultiplyAdd(a, b, c);
  } else {
    return c - a * b;
  }
}

// Evaluates the sum `x + Δx`.  If that sum has a dangerous rounding
// configuration (that is, the bits after the last mantissa bit of the sum are
// either 1000... or 0111..., then returns `NaN`.  Otherwise returns the sum.
// `x` is always positive.  `Δx` may be positive or negative.
inline double DetectDangerousRounding(double const x, double const Δx) {
  DoublePrecision<double> const sum = QuickTwoSum(x, Δx);
  double const& value = sum.value;
  double const& error = sum.error;
  __m128i const value_exponent_128i =
      _mm_castpd_si128(_mm_and_pd(masks::exponent_bits, _mm_set_sd(value)));
  __m128i const error_128i =
      _mm_castpd_si128(_mm_andnot_pd(masks::sign_bit, _mm_set_sd(error)));
  double const normalized_error = _mm_cvtsd_f64(
      _mm_castsi128_pd(_mm_sub_epi64(error_128i, value_exponent_128i)));
  // TODO(phl): Error analysis to refine the thresholds.  Also, can we avoid
  // negative numbers?
  OSACA_IF(normalized_error < -0x1.ffffp971 &&
           normalized_error > -0x1.00008p972) {
#if _DEBUG
    LOG(ERROR) << std::setprecision(25) << x << " " << std::hexfloat << value
               << " " << error << " " << normalized_error;
#endif
    return std::numeric_limits<double>::quiet_NaN();
  } else {
    return value;
  }
}

template<FMAPolicy fma_policy>
inline void Reduce(Argument const θ,
                   DoublePrecision<Argument>& θ_reduced,
                   std::int64_t& quadrant) {
  double const abs_θ = std::abs(θ);
  OSACA_IF(abs_θ < π / 4) {
    θ_reduced.value = θ;
    θ_reduced.error = 0;
    quadrant = 0;
    return;
  } OSACA_ELSE_IF(abs_θ <= π_over_2_threshold) {
    // We are not very sensitive to rounding errors in this expression, because
    // in the worst case it could cause the reduced angle to jump from the
    // vicinity of π / 4 to the vicinity of -π / 4 with appropriate adjustment
    // of the quadrant.
    __m128d const sign = _mm_and_pd(masks::sign_bit, _mm_set_sd(θ));
    double const n_shifted =
        FusedMultiplyAdd<fma_policy>(abs_θ, (2 / π), mantissa_reduce_shifter);
    double const n_double = _mm_cvtsd_f64(
        _mm_xor_pd(_mm_set_sd(n_shifted - mantissa_reduce_shifter), sign));
    std::int64_t const n = _mm_cvtsd_si64(_mm_set_sd(n_double));

    Argument const value =
        FusedNegatedMultiplyAdd<fma_policy>(n_double, π_over_2_high, θ);
    Argument const error = n_double * π_over_2_low;
    θ_reduced = QuickTwoDifference(value, error);
    // TODO(phl): Error analysis needed to find the right bounds.
    OSACA_IF(θ_reduced.value < -0x1.0p-30 || θ_reduced.value > 0x1.0p-30) {
      quadrant = n & 0b11;
      return;
    }
  }
  θ_reduced.value = 0;
  θ_reduced.error = std::numeric_limits<double>::quiet_NaN();
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
Value SinImplementation(DoublePrecision<Argument> const θ_reduced) {
  auto const& x = θ_reduced.value;
  auto const& e = θ_reduced.error;
  double const abs_x = std::abs(x);
  OSACA_IF(abs_x < sin_near_zero_cutoff) {
    double const x² = x * x;
    double const x³ = x² * x;
    double const x³_term = FusedMultiplyAdd<fma_policy>(
        x³, SinPolynomialNearZero<fma_policy>(x²), e);
    return DetectDangerousRounding(x, x³_term);
  } else {
    __m128d const sign = _mm_and_pd(masks::sign_bit, _mm_set_sd(x));
    double const e_abs = _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(e), sign));
    auto const i = AccurateTableIndex(abs_x);
    auto const& accurate_values = SinCosAccurateTable[i];
    double const& x₀ = accurate_values.x;
    double const& sin_x₀ = accurate_values.sin_x;
    double const& cos_x₀ = accurate_values.cos_x;
    // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
    // don't overlap and in the first interval `x` and `h` may be of the same
    // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
    // later in the computation.  Note that the terms in `e * h²` and higher are
    // *not* computed correctly below because they don't matter.
    double const h = abs_x - x₀;

    DoublePrecision<double> const sin_x₀_plus_h_cos_x₀ =
        TwoProductAdd<fma_policy>(cos_x₀, h, sin_x₀);
    double const h² = h * (h + 2 * e_abs);
    double const h³ = h² * h;
    double const polynomial_term =
        FusedMultiplyAdd<fma_policy>(
            cos_x₀,
            h³ * SinPolynomial<fma_policy>(h²),
            sin_x₀ * h² * CosPolynomial<fma_policy>(h²)) +
        FusedMultiplyAdd<fma_policy>(cos_x₀, e_abs, sin_x₀_plus_h_cos_x₀.error);
    return _mm_cvtsd_f64(
        _mm_xor_pd(_mm_set_sd(DetectDangerousRounding(
                       sin_x₀_plus_h_cos_x₀.value, polynomial_term)),
                   sign));
  }
}

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value CosImplementation(DoublePrecision<Argument> const θ_reduced) {
  auto const& x = θ_reduced.value;
  auto const& e = θ_reduced.error;
  double const abs_x = std::abs(x);
  __m128d const sign = _mm_and_pd(masks::sign_bit, _mm_set_sd(x));
  double const e_abs = _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(e), sign));
  auto const i = AccurateTableIndex(abs_x);
  auto const& accurate_values = SinCosAccurateTable[i];
  double const& x₀ = accurate_values.x;
  double const& sin_x₀ = accurate_values.sin_x;
  double const& cos_x₀ = accurate_values.cos_x;
  // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
  // don't overlap and in the first interval `x` and `h` may be of the same
  // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
  // later in the computation.  Note that the terms in `e * h²` and higher are
  // *not* computed correctly below because they don't matter.
  double const h = abs_x - x₀;

  DoublePrecision<double> const cos_x₀_minus_h_sin_x₀ =
      TwoProductNegatedAdd<fma_policy>(sin_x₀, h, cos_x₀);
  double const h² = h * (h + 2 * e_abs);
  double const h³ = h² * h;
  double const polynomial_term =
      FusedNegatedMultiplyAdd<fma_policy>(
          sin_x₀,
          h³ * SinPolynomial<fma_policy>(h²),
          cos_x₀ * h² * CosPolynomial<fma_policy>(h²)) +
      FusedNegatedMultiplyAdd<fma_policy>(
          sin_x₀, e_abs, cos_x₀_minus_h_sin_x₀.error);
  return DetectDangerousRounding(cos_x₀_minus_h_sin_x₀.value, polynomial_term);
}

Value __cdecl Sin(Argument θ) {
  OSACA_FUNCTION_BEGIN(θ);
  DoublePrecision<Argument> θ_reduced;
  std::int64_t quadrant;
  double value;
  OSACA_IF(UseHardwareFMA) {
    Reduce<FMAPolicy::Force>(θ, θ_reduced, quadrant);
    OSACA_IF(quadrant & 0b1) {
      value = CosImplementation<FMAPolicy::Force>(θ_reduced);
    } else {
      value = SinImplementation<FMAPolicy::Force>(θ_reduced);
    }
  } else {
    Reduce<FMAPolicy::Disallow>(θ, θ_reduced, quadrant);
    OSACA_IF(quadrant & 0b1) {
      value = CosImplementation<FMAPolicy::Disallow>(θ_reduced);
    } else {
      value = SinImplementation<FMAPolicy::Disallow>(θ_reduced);
    }
  }
  OSACA_IF(value != value) {
    OSACA_RETURN(cr_sin(θ));
  } OSACA_ELSE_IF(quadrant & 0b10) {
    OSACA_RETURN(-value);
  } else {
    OSACA_RETURN(value);
  }
}

Value __cdecl Cos(Argument θ) {
  OSACA_FUNCTION_BEGIN(θ);
  DoublePrecision<Argument> θ_reduced;
  std::int64_t quadrant;
  double value;
  OSACA_IF(UseHardwareFMA) {
    Reduce<FMAPolicy::Force>(θ, θ_reduced, quadrant);
    OSACA_IF(quadrant & 0b1) {
      value = SinImplementation<FMAPolicy::Force>(θ_reduced);
    } else {
      value = CosImplementation<FMAPolicy::Force>(θ_reduced);
    }
  } else {
    Reduce<FMAPolicy::Disallow>(θ, θ_reduced, quadrant);
    OSACA_IF(quadrant & 0b1) {
      value = SinImplementation<FMAPolicy::Disallow>(θ_reduced);
    } else {
      value = CosImplementation<FMAPolicy::Disallow>(θ_reduced);
    }
  }
  OSACA_IF(value != value) {
    OSACA_RETURN(cr_cos(θ));
  } OSACA_ELSE_IF(quadrant == 1 || quadrant == 2) {
    OSACA_RETURN(-value);
  } else {
    OSACA_RETURN(value);
  }
}

}  // namespace internal
}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
