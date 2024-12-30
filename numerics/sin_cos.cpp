#pragma once

#include "numerics/sin_cos.hpp"

#include <pmmintrin.h>

#include <limits>

#include "core-math/cos.h"
#include "core-math/sin.h"
#include "numerics/accurate_tables.mathematica.h"
#include "numerics/double_precision.hpp"
#include "numerics/fma.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"

#if PRINCIPIA_USE_OSACA_SIN || PRINCIPIA_USE_OSACA_COS
#include "intel/iacaMarks.h"
static bool const iaca_loop_terminator = false;
#define IACA_FUNCTION_DOUBLE(arg) \
  double iaca_result = arg;       \
  IACA_VC64_START;                \
  IACA_LOOP :                    \
    arg = iaca_result
#define IACA_RETURN(result) IACA_RETURN_AT(__LINE__, result)
#define IACA_RETURN_AT(line, result) IACA_RETURN_AT2(line, result)
#define IACA_RETURN_AT2(line, result)  \
  if (iaca_loop_terminator) {          \
    goto IACA_END_LOOP_##line;          \
  }                                    \
  iaca_result = (result);              \
  goto IACA_LOOP;                      \
  IACA_VC64_END; \
  IACA_END_LOOP_##line : return iaca_result
#define IACA_FUNCTION_END
  #define IACA_IF(condition)                                           \
  if constexpr (volatile bool IACA_computed_condition = (condition); \
                [] { UNDER_IACA_HYPOTHESES(return (condition)); }())

#define UNDER_IACA_HYPOTHESES(statement)                                       \
  do {                                                                         \
    constexpr bool UseHardwareFMA = true;                                      \
    constexpr double θ = 0.1;                                                  \
    /* From argument reduction. */                                             \
    constexpr double n_double = θ * (2 / π);                                   \
    constexpr double reduction_value = θ - n_double * π_over_2_high;           \
    constexpr double reduction_error = n_double * π_over_2_low;                \
    /* Used to determine whether a better argument reduction is needed. */     \
    constexpr DoublePrecision<double> θ_reduced =                              \
        QuickTwoDifference(reduction_value, reduction_error);                  \
    /* Used in Sin to detect the near-0 case. */                               \
    constexpr double abs_x =                                                   \
        θ_reduced.value > 0 ? θ_reduced.value : -θ_reduced.value;              \
    /* Used throughout the top-level functions. */                             \
    constexpr std::int64_t quadrant =                                          \
        static_cast<std::int64_t>(n_double) & 0b11;                            \
    /* Used in DetectDangerousRounding. */                                     \
    constexpr double normalized_error = 0;                                     \
    /* Not NaN is the only part that matters; used at the end of the top-level \
     * functions to determine whether to call the slow path. */                \
    constexpr double value = 1;                                                \
    { statement; }                                                             \
  } while (false)

#else
#define IACA_FUNCTION_DOUBLE(arg) (double const arg)
#define IACA_RETURN(result) return result
#define IACA_IF (condition) if (condition)
#define UNDER_IACA_HYPOTHESES(statement)
#endif


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

constexpr std::int64_t table_spacing_bits = 9;
constexpr Argument table_spacing_reciprocal = 1 << table_spacing_bits;
constexpr Argument table_spacing = 1.0 / table_spacing_reciprocal;
constexpr Argument sin_near_zero_cutoff = table_spacing / 2.0;

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
  IACA_IF ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
          (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using quantities::_elementary_functions::FusedMultiplyAdd;
    return FusedMultiplyAdd(a, b, c);
  } else {
    return a * b + c;
  }
}

template<FMAPolicy fma_policy>
double FusedNegatedMultiplyAdd(double const a, double const b, double const c) {
  IACA_IF ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
          (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using quantities::_elementary_functions::FusedNegatedMultiplyAdd;
    return FusedNegatedMultiplyAdd(a, b, c);
  } else {
    return c - a * b;
  }
}

// Evaluates the sum `x + Δx`.  If that sum has a dangerous rounding
// configuration (that is, the bits after the last mantissa bit of the sum are
// either 1000... or 0111..., then returns `NaN`.  Otherwise returns the sum.
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
  IACA_IF (normalized_error < -0x1.ffffp971 &&
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

inline void Reduce(Argument const θ,
                   DoublePrecision<Argument>& θ_reduced,
                   std::int64_t& quadrant) {
  IACA_IF (θ < π / 4 && θ > -π / 4) {
    θ_reduced.value = θ;
    θ_reduced.error = 0;
    quadrant = 0;
    return;
  } else IACA_IF (θ <= π_over_2_threshold && θ >= -π_over_2_threshold) {
    // We are not very sensitive to rounding errors in this expression, because
    // in the worst case it could cause the reduced angle to jump from the
    // vicinity of π / 4 to the vicinity of -π / 4 with appropriate adjustment
    // of the quadrant.
    __m128d const n_128d = _mm_round_sd(
        _mm_setzero_pd(), _mm_set_sd(θ * (2 / π)), _MM_FROUND_RINT);
    double const n_double = _mm_cvtsd_f64(n_128d);
    std::int64_t const n = _mm_cvtsd_si64(n_128d);
    Argument const value = θ - n_double * π_over_2_high;
    Argument const error = n_double * π_over_2_low;
    θ_reduced = QuickTwoDifference(value, error);
    // TODO(phl): Error analysis needed to find the right bounds.
    IACA_IF (θ_reduced.value < -0x1.0p-30 || θ_reduced.value > 0x1.0p-30) {
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
  IACA_IF (abs_x < sin_near_zero_cutoff) {
    double const x² = x * x;
    double const x³ = x² * x;
    double const x³_term = FusedMultiplyAdd<fma_policy>(
                   x³, SinPolynomialNearZero<fma_policy>(x²), e);
    return DetectDangerousRounding(x, x³_term);
  } else {
    __m128d const sign = _mm_and_pd(masks::sign_bit, _mm_set_sd(x));
    auto const i = AccurateTableIndex(abs_x);
    auto const& accurate_values = SinCosAccurateTable[i];
    double const x₀ =
        _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.x), sign));
    double const sin_x₀ =
        _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.sin_x), sign));
    double const& cos_x₀ = accurate_values.cos_x;
    // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
    // don't overlap and in the first interval `x` and `h` may be of the same
    // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
  // later in the computation.  Note that the terms in `e * h²` and higher are
  // *not* computed correctly below because they don't matter.
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
    return DetectDangerousRounding(sin_x₀_plus_h_cos_x₀.value, polynomial_term);
  }
}

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value CosImplementation(DoublePrecision<Argument> const θ_reduced) {
  auto const& x = θ_reduced.value;
  auto const& e = θ_reduced.error;
  double const abs_x = std::abs(x);
  __m128d const sign = _mm_and_pd(masks::sign_bit, _mm_set_sd(x));
  auto const i = AccurateTableIndex(abs_x);
  auto const& accurate_values = SinCosAccurateTable[i];
  double const x₀ =
      _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.x), sign));
  double const sin_x₀ =
      _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(accurate_values.sin_x), sign));
  double const& cos_x₀ = accurate_values.cos_x;
  // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
  // don't overlap and in the first interval `x` and `h` may be of the same
  // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
  // later in the computation.  Note that the terms in `e * h²` and higher are
  // *not* computed correctly below because they don't matter.
  double const h = x - x₀;

  DoublePrecision<double> const cos_x₀_minus_h_sin_x₀ =
      TwoProductNegatedAdd<fma_policy>(sin_x₀, h, cos_x₀);
  double const h² = h * (h + 2 * e);
  double const h³ = h² * h;
  double const polynomial_term =
      (FusedNegatedMultiplyAdd<fma_policy>(
           sin_x₀,
           FusedMultiplyAdd<fma_policy>(h³, SinPolynomial<fma_policy>(h²), e),
           cos_x₀ * h² * CosPolynomial<fma_policy>(h²)) +
       cos_x₀_minus_h_sin_x₀.error);
  return DetectDangerousRounding(cos_x₀_minus_h_sin_x₀.value, polynomial_term);
}

#if PRINCIPIA_INLINE_SIN_COS
FORCE_INLINE(inline)
#endif
Value __cdecl Sin(Argument const θ) {
  DoublePrecision<Argument> θ_reduced;
  std::int64_t quadrant;
  Reduce(θ, θ_reduced, quadrant);
  double value;
  IACA_IF (UseHardwareFMA) {
    IACA_IF (quadrant & 0b1) {
      value = CosImplementation<FMAPolicy::Force>(θ_reduced);
    } else {
#if PRINCIPIA_USE_OSACA_SIN
      IACA_VC64_START;
#endif
      value = SinImplementation<FMAPolicy::Force>(θ_reduced);
#if PRINCIPIA_USE_OSACA_SIN
      IACA_VC64_END;
#endif
    }
  } else {
    IACA_IF (quadrant & 0b1) {
      value = CosImplementation<FMAPolicy::Disallow>(θ_reduced);
    } else {
      value = SinImplementation<FMAPolicy::Disallow>(θ_reduced);
    }
  }
  IACA_IF (value != value) {
    return cr_sin(θ);
  }
  IACA_IF (quadrant & 0b10) {
    return -value;
  } else {
    return value;
  }
}

namespace IACA_assumptions {
constexpr std::int64_t quadrant = 1;
}  // namespace IACA_assumptions

#if PRINCIPIA_INLINE_SIN_COS
FORCE_INLINE(inline)
#endif
Value __cdecl Cos(Argument θ) {
  IACA_FUNCTION_DOUBLE(θ);
  DoublePrecision<Argument> θ_reduced;
  std::int64_t quadrant;
  Reduce(θ, θ_reduced, quadrant);
  double value;
  IACA_IF (UseHardwareFMA) {
    IACA_IF (quadrant & 0b1) {
      value = SinImplementation<FMAPolicy::Force>(θ_reduced);
    } else {
      value = CosImplementation<FMAPolicy::Force>(θ_reduced);
    }
  } else {
    IACA_IF (quadrant & 0b1) {
      value = SinImplementation<FMAPolicy::Disallow>(θ_reduced);
    } else {
      value = CosImplementation<FMAPolicy::Disallow>(θ_reduced);
    }
  }
  IACA_IF (value != value) {
    IACA_RETURN(cr_cos(θ));
  } else IACA_IF (quadrant == 1 || quadrant == 2) {
    IACA_RETURN(-value);
  } else {
    IACA_RETURN(value);
  }
  IACA_FUNCTION_END;
}

}  // namespace internal
}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
