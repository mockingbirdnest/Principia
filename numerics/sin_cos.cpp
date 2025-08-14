#pragma once

#include "numerics/sin_cos.hpp"

#include <pmmintrin.h>

#include <limits>

#include "core-math/cos.h"
#include "core-math/sin.h"
#include "numerics/accurate_tables.mathematica.h"
#include "numerics/double_precision.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/fma.hpp"
#include "numerics/osaca.hpp"  // 🧙 For OSACA_*.
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/m128d.hpp"

// The algorithms in this file are documented in `Sin Cos.pdf`.  To the extent
// possible, the code follows the notation of that document.
namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

using namespace principia::numerics::_accurate_tables;
using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_fma;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::quantities::_m128d;

#define OSACA_ANALYSED_FUNCTION Cos
#define OSACA_ANALYSED_FUNCTION_NAMESPACE
#define OSACA_ANALYSED_FUNCTION_TEMPLATE_PARAMETERS <FMAPolicy::Force>
#define UNDER_OSACA_HYPOTHESES(expression)                                   \
  [&] {                                                                      \
    constexpr bool UseHardwareFMA = true;                                    \
    constexpr double θ = 3;                                                  \
    /* From argument reduction. */                                           \
    constexpr double abs_θ = θ > 0 ? θ : -θ;                                 \
    constexpr std::int64_t n = static_cast<std::int64_t>(θ * (2 / π) + 0.5); \
    constexpr double reduction_value = θ - n * C₁;                           \
    constexpr double reduction_error = n * δC₁;                              \
    /* Used to determine whether a better argument reduction is needed. */   \
    constexpr DoublePrecision<double> θ_reduced =                            \
        TwoDifference(reduction_value, reduction_error);                     \
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

constexpr std::int64_t κ₁ = 8;
constexpr std::int64_t κʹ₁ = 5;
constexpr std::int64_t κ₂ = 18;
constexpr std::int64_t κʹ₂ = 14;
constexpr std::int64_t κʺ₂ = 15;
constexpr std::int64_t κ₃ = 18;
constexpr Argument two_term_θ_threshold = π / 2 * (1LL << κ₁);
constexpr Argument three_term_θ_threshold = π / 2 * (1LL << κ₂);
constexpr Argument C₁ = 0x1.921F'B544'42D0'0p0;
constexpr Argument δC₁ = 0x1.8469'898C'C517'0p-48;
constexpr Argument C₂ = 0x1.921F'B544'4000'0p0;
constexpr Argument Cʹ₂ = 0x1.68C2'34C4'C000'0p-39;
constexpr Argument δC₂ = 0x1.98A2'E037'0734'5p-77;
constexpr Argument two_term_θ_reduced_threshold =
    1.0 / (1LL << (-(κ₁ + κʹ₁ + κ₃ - std::numeric_limits<double>::digits + 2)));
constexpr Argument three_term_θ_reduced_threshold =
    (1.0 / (1LL << (-(κ₃ - std::numeric_limits<double>::digits)))) *
    ((1LL << (-(κ₂ + κʹ₂ + κʺ₂ - std::numeric_limits<double>::digits + 2))) +
     4);

constexpr Argument mantissa_reduce_shifter =
    1LL << (std::numeric_limits<double>::digits - 1);

constexpr double sin_near_zero_e = 0x1.0000'32D7'5E64'Cp0;
constexpr double sin_e = 0x1.0000'A03C'34D3'9p0;
constexpr double cos_e = 0x1.0000'A07E'1980'8p0;

template<FMAPolicy fma_policy>
using Polynomial1 = HornerEvaluator<M128D, M128D, 1, fma_policy>;

// Pointers used for indirect calls, set by `StaticInitialization`.
Value (__cdecl *cos)(Argument θ) = nullptr;
Value (__cdecl *sin)(Argument θ) = nullptr;

// Forward declarations needed by the OSACA macros.
template<FMAPolicy fma_policy>
Value __cdecl Sin(Argument θ);
template<FMAPolicy fma_policy>
Value __cdecl Cos(Argument θ);

namespace masks {
M128D const sign_bit(0x8000'0000'0000'0000ull);
M128D const exponent_bits(0x7ff0'0000'0000'0000ull);
M128D const mantissa_bits(0x000f'ffff'ffff'ffffull);
M128D const mantissa_index_bits(0x0000'0000'0000'01ffull);
}  // namespace masks

inline std::int64_t AccurateTableIndex(M128D const abs_x) {
  // This function computes the index in the accurate table:
  // 1. A suitable (large) power of 2 is added to the argument so that the last
  //    bit of the mantissa of the result corresponds to units of 1/512 of the
  //    argument.  As part of this addition, an interval of radius 1/1024 around
  //    an integral multiple of 1/512 is correctly rounded to that integral
  //    multiple.
  // 2. An `and` operation is used to only retain the last 9 bits of the
  //    mantissa.
  // 3. The result is interpreted as an integer and returned as the index.
  return static_cast<std::int64_t>(
      masks::mantissa_index_bits &
      abs_x + M128D(1LL << (std::numeric_limits<double>::digits -
                            table_spacing_bits - 1)));
}

template<FMAPolicy fma_policy>
M128D FusedMultiplyAdd(M128D const a, M128D const b, M128D const c) {
  static_assert(fma_policy != FMAPolicy::Auto);
  if constexpr (fma_policy == FMAPolicy::Force) {
    return FusedMultiplyAdd(a, b, c);
  } else {
    return a * b + c;
  }
}

template<FMAPolicy fma_policy>
M128D FusedNegatedMultiplyAdd(M128D const a, M128D const b, M128D const c) {
  static_assert(fma_policy != FMAPolicy::Auto);
  if constexpr (fma_policy == FMAPolicy::Force) {
    return FusedNegatedMultiplyAdd(a, b, c);
  } else {
    return c - a * b;
  }
}

template<FMAPolicy fma_policy>
double FusedMultiplyAdd(double const a, double const b, double const c) {
  return static_cast<double>(
      FusedMultiplyAdd<fma_policy>(M128D(a), M128D(b), M128D(c)));
}

template<FMAPolicy fma_policy>
double FusedNegatedMultiplyAdd(double const a, double const b, double const c) {
  return static_cast<double>(
      FusedNegatedMultiplyAdd<fma_policy>(M128D(a), M128D(b), M128D(c)));
}

// Evaluates the sum `x + Δx` and performs the rounding test using the technique
// described in [Mul+10], section 11.6.3.  If rounding is safe, returns the sum;
// otherwise, returns `NaN`.  `x` is always positive.  `Δx` may be positive or
// negative.
template<FMAPolicy fma_policy, double e>
double DetectDangerousRounding(M128D const x, M128D const Δx) {
  // We don't check that `Δx` is not NaN because that's how we trigger fallback
  // to the slow path.
  DCHECK(x == x);
  DoublePrecision<M128D> const sum = QuickTwoSum(x, Δx);
  auto const& value = sum.value;
  auto const& error = sum.error;
  OSACA_IF(value == FusedMultiplyAdd<fma_policy>(error, e, value)) {
    return static_cast<double>(value);
  } else {
#if _DEBUG
    LOG_IF(ERROR, value == value && error == error)
        << std::setprecision(25) << x << " " << std::hexfloat << value << " "
        << error << " " << e;
#endif
    return std::numeric_limits<double>::quiet_NaN();
  }
}

template<FMAPolicy fma_policy, bool preserve_sign>
FORCE_INLINE(inline)
void Reduce(Argument const θ,
            DoublePrecision<Argument>& θ_reduced,
            std::int64_t& quadrant) {
  double const abs_θ = std::abs(θ);
  OSACA_IF(abs_θ < π / 4) {
    θ_reduced.value = θ;
    θ_reduced.error = 0;
    quadrant = 0;
    return;
  } OSACA_ELSE_IF(abs_θ <= two_term_θ_threshold) {
    // We are not very sensitive to rounding errors in this expression, because
    // in the worst case it could cause the reduced angle to jump from the
    // vicinity of π / 4 to the vicinity of -π / 4 with appropriate adjustment
    // of the quadrant.
    M128D const sign = Sign(θ);
    double n_double =
        FusedMultiplyAdd<fma_policy>(abs_θ, (2 / π), mantissa_reduce_shifter) -
        mantissa_reduce_shifter;

    // Don't move the computation of `n` after the if, it generates some extra
    // moves.
    Argument y;
    std::int64_t n;
    if constexpr (preserve_sign) {
      n_double = static_cast<double>(n_double ^ sign);
      n = _mm_cvtsd_si64(_mm_set_sd(n_double));
      y = FusedNegatedMultiplyAdd<fma_policy>(n_double, C₁, θ);
    } else {
      n = _mm_cvtsd_si64(_mm_set_sd(n_double));
      y = FusedNegatedMultiplyAdd<fma_policy>(n_double, C₁, abs_θ);
    }

    Argument const δy = n_double * δC₁;
    θ_reduced = TwoDifference(y, δy);
    OSACA_IF(θ_reduced.value <= -two_term_θ_reduced_threshold ||
             θ_reduced.value >= two_term_θ_reduced_threshold) {
      quadrant = n & 0b11;
      return;
    }
  } OSACA_ELSE_IF(abs_θ <= three_term_θ_threshold) {
    // Same code as above.
    M128D const sign = Sign(θ);
    double n_double =
        FusedMultiplyAdd<fma_policy>(abs_θ, (2 / π), mantissa_reduce_shifter) -
        mantissa_reduce_shifter;

    Argument y;
    std::int64_t n;
    if constexpr (preserve_sign) {
      n_double = static_cast<double>(n_double ^ sign);
      n = _mm_cvtsd_si64(_mm_set_sd(n_double));
      y = FusedNegatedMultiplyAdd<fma_policy>(n_double, C₂, θ);
    } else {
      n = _mm_cvtsd_si64(_mm_set_sd(n_double));
      y = FusedNegatedMultiplyAdd<fma_policy>(n_double, C₂, abs_θ);
    }

    Argument const yʹ = n_double * Cʹ₂;
    Argument const δy = n_double * δC₂;
    auto const z = QuickTwoSum(yʹ, δy);
    θ_reduced = y - z;
    OSACA_IF(θ_reduced.value <= -three_term_θ_reduced_threshold ||
             θ_reduced.value >= three_term_θ_reduced_threshold) {
      quadrant = n & 0b11;
      return;
    }
  }
  θ_reduced.value = 0;
  θ_reduced.error = std::numeric_limits<double>::quiet_NaN();
}

template<FMAPolicy fma_policy>
M128D SinPolynomial(M128D const x) {
  // Absolute error of the exact polynomial better than 85.7 bits over an
  // interval slightly larger than [-1/1024, 1/1024].
  return Polynomial1<fma_policy>::Evaluate(
      {-0x1.5555'5555'5555'5p-3, 0x1.1111'10B1'75B5'Fp-7}, x);
}

template<FMAPolicy fma_policy>
M128D SinPolynomialNearZero(M128D const x) {
  // Relative error of the exact polynomial better than 75.5 bits over
  // [-1/1024, 1/1024].
  return Polynomial1<fma_policy>::Evaluate(
      {-0x1.5555'5555'5555'5p-3, 0x1.1111'10B4'0E88'Ap-7}, x);
}

template<FMAPolicy fma_policy>
M128D CosPolynomial(M128D const x) {
  // Absolute error of the exact polynomial better than 72.6 bits over an
  // interval slightly larger than [-1/1024, 1/1024].
  return Polynomial1<fma_policy>::Evaluate(
      {-0x1.0000'0000'0000'0p-1, 0x1.5555'54B1'22F2'9p-5}, x);
}

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value SinImplementation(DoublePrecision<Argument> const θ_reduced) {
  M128D const x = θ_reduced.value;
  M128D const e = θ_reduced.error;
  auto const abs_x = Abs(x);
  OSACA_IF(abs_x < sin_near_zero_cutoff) {
    auto const x² = x * x;
    auto const x³ = x² * x;
    auto const x³_term = FusedMultiplyAdd<fma_policy>(
        x³, SinPolynomialNearZero<fma_policy>(x²), e);
    // Relative error of the result better than 72.8 bits.
    return DetectDangerousRounding<fma_policy, sin_near_zero_e>(x, x³_term);
  } else {
    auto const sign = Sign(x);
    auto const e_abs = e ^ sign;
    auto const i = AccurateTableIndex(abs_x);
    auto const& accurate_values = SinCosAccurateTable[i];
    M128D const x₀ = accurate_values.x;
    M128D const sin_x₀ = accurate_values.sin_x;
    M128D const cos_x₀ = accurate_values.cos_x;
    // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
    // don't overlap and in the first interval `x` and `h` may be of the same
    // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
    // later in the computation.  Note that the terms in `e * h²` and higher are
    // *not* computed correctly below because they don't matter.
    auto const h = abs_x - x₀;

    DoublePrecision<M128D> const sin_x₀_plus_h_cos_x₀ =
        TwoProductAdd<fma_policy>(cos_x₀, h, sin_x₀);
    auto const h² = h * (h + 2 * e_abs);
    auto const h³ = h² * h;
    auto const polynomial_term =
        FusedMultiplyAdd<fma_policy>(
            cos_x₀,
            h³ * SinPolynomial<fma_policy>(h²),
            sin_x₀ * h² * CosPolynomial<fma_policy>(h²)) +
        FusedMultiplyAdd<fma_policy>(cos_x₀, e_abs, sin_x₀_plus_h_cos_x₀.error);
    return static_cast<double>(
        sign ^ DetectDangerousRounding<fma_policy, sin_e>(
                   sin_x₀_plus_h_cos_x₀.value, polynomial_term));
  }
}

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value CosImplementation(DoublePrecision<Argument> const θ_reduced) {
  M128D const x = θ_reduced.value;
  M128D const e = θ_reduced.error;
  auto const abs_x = Abs(x);
  auto const sign = Sign(x);
  auto const e_abs = e ^ sign;
  auto const i = AccurateTableIndex(abs_x);
  auto const& accurate_values = SinCosAccurateTable[i];
  M128D const x₀ = accurate_values.x;
  M128D const sin_x₀ = accurate_values.sin_x;
  M128D const cos_x₀ = accurate_values.cos_x;
  // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
  // don't overlap and in the first interval `x` and `h` may be of the same
  // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
  // later in the computation.  Note that the terms in `e * h²` and higher are
  // *not* computed correctly below because they don't matter.
  auto const h = abs_x - x₀;

  DoublePrecision<M128D> const cos_x₀_minus_h_sin_x₀ =
      TwoProductNegatedAdd<fma_policy>(sin_x₀, h, cos_x₀);
  auto const h² = h * (h + 2 * e_abs);
  auto const h³ = h² * h;
  auto const polynomial_term =
      FusedNegatedMultiplyAdd<fma_policy>(
          sin_x₀,
          h³ * SinPolynomial<fma_policy>(h²),
          cos_x₀ * h² * CosPolynomial<fma_policy>(h²)) +
      FusedNegatedMultiplyAdd<fma_policy>(
          sin_x₀, e_abs, cos_x₀_minus_h_sin_x₀.error);
  return DetectDangerousRounding<fma_policy, cos_e>(cos_x₀_minus_h_sin_x₀.value,
                                                    polynomial_term);
}

template<FMAPolicy fma_policy>
Value __cdecl Sin(Argument θ) {
  OSACA_FUNCTION_BEGIN(θ, <fma_policy>);
  DoublePrecision<Argument> θ_reduced;
  std::int64_t quadrant;
  double value;
  Reduce<fma_policy, /*preserve_sign=*/true>(θ, θ_reduced, quadrant);
  OSACA_IF(quadrant & 0b1) {
    value = CosImplementation<fma_policy>(θ_reduced);
  } else {
    value = SinImplementation<fma_policy>(θ_reduced);
  }
  OSACA_IF(value != value) {
    OSACA_RETURN(cr_sin(θ));
  } OSACA_ELSE_IF(quadrant & 0b10) {
    OSACA_RETURN(-value);
  } else {
    OSACA_RETURN(value);
  }
}

template<FMAPolicy fma_policy>
Value __cdecl Cos(Argument θ) {
  OSACA_FUNCTION_BEGIN(θ, <fma_policy>);
  DoublePrecision<Argument> θ_reduced;
  std::int64_t quadrant;
  double value;
  Reduce<fma_policy, /*preserve_sign=*/false>(θ, θ_reduced, quadrant);
  OSACA_IF(quadrant & 0b1) {
    value = SinImplementation<fma_policy>(θ_reduced);
  } else {
    value = CosImplementation<fma_policy>(θ_reduced);
  }
  OSACA_IF(value != value) {
    OSACA_RETURN(cr_cos(θ));
  } OSACA_ELSE_IF(quadrant == 1 || quadrant == 2) {
    OSACA_RETURN(-value);
  } else {
    OSACA_RETURN(value);
  }
}

void StaticInitialization() {
  if (UseHardwareFMA) {
    cos = &Cos<FMAPolicy::Force>;
    sin = &Sin<FMAPolicy::Force>;
  } else {
    cos = &Cos<FMAPolicy::Disallow>;
    sin = &Sin<FMAPolicy::Disallow>;
  }
}

Value __cdecl Sin(Argument const θ) {
  return sin(θ);
}

Value __cdecl Cos(Argument const θ) {
  return cos(θ);
}

}  // namespace internal
}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
