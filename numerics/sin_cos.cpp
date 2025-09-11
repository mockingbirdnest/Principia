#include "numerics/sin_cos.hpp"

#include <pmmintrin.h>

#include <limits>
#include <utility>

#include "base/tags.hpp"
#include "core-math/cos.h"
#include "core-math/sin.h"
#include "numerics/accurate_tables.mathematica.h"
#include "numerics/double_precision.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/osaca.hpp"  // 🧙 For OSACA_*.
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/m128d.hpp"

// The algorithms in this file are documented in `Sin Cos.pdf`.  To the extent
// possible, the code follows the notation of that document.
namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::numerics::_accurate_tables;
using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_m128d;

#define OSACA_ANALYSED_FUNCTION Cos
#define OSACA_ANALYSED_FUNCTION_NAMESPACE
#if PRINCIPIA_COMPILER_MSVC
#define OSACA_ANALYSED_FUNCTION_TEMPLATE_PARAMETERS <FMAPresence::Present>
#else
#define OSACA_ANALYSED_FUNCTION_TEMPLATE_PARAMETERS <FMAPresence::Absent>
#endif
#define UNDER_OSACA_HYPOTHESES(expression)                                   \
  [&] {                                                                      \
    constexpr bool CanUseHardwareFMA = true;                                 \
    constexpr double x = 3;                                                  \
    /* From argument reduction. */                                           \
    constexpr double abs_x = x > 0 ? x : -x;                                 \
    constexpr std::int64_t n = static_cast<std::int64_t>(x * (2 / π) + 0.5); \
    constexpr double reduction_value = x - n * C₁;                           \
    constexpr double reduction_error = n * δC₁;                              \
    /* Used to determine whether a better argument reduction is needed. */   \
    constexpr DoublePrecision<double> x_reduced =                            \
        TwoDifference(reduction_value, reduction_error);                     \
    constexpr double abs_x_reduced_value =                                   \
        x_reduced.value > 0 ? x_reduced.value : -x_reduced.value;            \
    /* Used in Sin to detect the near-0 case. */                             \
    constexpr double abs_x̃ =                                                 \
        x_reduced.value > 0 ? x_reduced.value : -x_reduced.value;            \
    /* Used throughout the top-level functions. */                           \
    constexpr std::int64_t quadrant = n & 0b11;                              \
    /* Not NaN is the only part that matters; used at the end of the    */   \
    /* top-level functions to determine whether to call the slow path.  */   \
    constexpr double value = 1;                                              \
    constexpr double muller_test_expression = value;                         \
    constexpr SC<double> values = {.sin = 0, .cos = 1};                      \
    return expression;                                                       \
  }()

using Argument = M128D;
using Value = M128D;

constexpr std::int64_t table_spacing_bits = 9;
constexpr double table_spacing_reciprocal = 1 << table_spacing_bits;
constexpr double table_spacing = 1.0 / table_spacing_reciprocal;
constexpr double sin_near_zero_cutoff = table_spacing / 2.0;

constexpr std::int64_t κ₁ = 8;
constexpr std::int64_t κʹ₁ = 5;
constexpr std::int64_t κ₂ = 18;
constexpr std::int64_t κʹ₂ = 14;
constexpr std::int64_t κʺ₂ = 15;
constexpr std::int64_t κ₃ = 18;

// These constants must be `constexpr` (and therefore `double`) to be used in
// the `OSACA_` macros.
constexpr double two_term_x_threshold = π / 2 * (1LL << κ₁);
constexpr double three_term_x_threshold = π / 2 * (1LL << κ₂);
constexpr double C₁ = 0x1.921F'B544'42D0'0p0;
constexpr double δC₁ = 0x1.8469'898C'C517'0p-48;
constexpr double C₂ = 0x1.921F'B544'4000'0p0;
constexpr double Cʹ₂ = 0x1.68C2'34C4'C000'0p-39;
constexpr double δC₂ = 0x1.98A2'E037'0734'5p-77;
constexpr double two_term_x_reduced_threshold =
    1.0 / (1LL << (-(κ₁ + κʹ₁ + κ₃ - std::numeric_limits<double>::digits + 2)));
constexpr double three_term_x_reduced_threshold =
    (1.0 / (1LL << (-(κ₃ - std::numeric_limits<double>::digits)))) *
    ((1LL << (-(κ₂ + κʹ₂ + κʺ₂ - std::numeric_limits<double>::digits + 2))) +
     4);

constexpr double sin_near_zero_e = 0x1.0000'B2D8'01D8'0p0;  // 2^-70.517.
constexpr double sin_e = 0x1.0002'6ABE'272D'9p0;  // 2^-68.726.
constexpr double cos_e = 0x1.0001'B839'335E'1p0;  // 2^-69.217.

SlowPathCallback slow_path_sin_callback = nullptr;
SlowPathCallback slow_path_cos_callback = nullptr;

// Forward declarations needed by the OSACA macros.
template<FMAPresence fma_presence>
double __cdecl Sin(double x);
template<FMAPresence fma_presence>
double __cdecl Cos(double x);

namespace m128d {

M128D const quiet_NaN(std::numeric_limits<double>::quiet_NaN());
M128D const zero(0.0);

// Argument reduction.
M128D const mantissa_reduce_shifter(
    static_cast<double>(1LL << (std::numeric_limits<double>::digits - 1)));
M128D const two_over_π(2.0 / π);
M128D const C₁(internal::C₁);
M128D const δC₁(internal::δC₁);
M128D const C₂(internal::C₂);
M128D const Cʹ₂(internal::Cʹ₂);
M128D const δC₂(internal::δC₂);

// Accurate table index.
M128D const mantissa_index_bits = M128D::MakeFromBits(0x0000'0000'0000'01ffull);
M128D const accurate_table_index_addend(static_cast<double>(
    1LL << (std::numeric_limits<double>::digits - table_spacing_bits - 1)));

// Polynomials.
M128D const sin_0(-0x1.5555'5555'5555'5p-3);
M128D const sin_1(0x1.1111'10A8'20AE'Cp-7);
M128D const sin_near_zero_0(-0x1.5555'5555'5555'5p-3);
M128D const sin_near_zero_1(0x1.1111'10B4'0E88'Ap-7);
M128D const cos_0(-0x1.FFFF'FFFF'FFFF'Dp-2);
M128D const cos_1(0x1.5555'549D'B0A9'5p-5);

}  // namespace m128d

template<FMAPresence fma_presence>
M128D MaybeFusedMultiplyAdd(M128D const a, M128D const b, M128D const c) {
  static_assert(fma_presence != FMAPresence::Unknown);
  if constexpr (fma_presence == FMAPresence::Present) {
    return FusedMultiplyAdd(a, b, c);
  } else {
    return a * b + c;
  }
}

template<FMAPresence fma_presence>
M128D MaybeFusedNegatedMultiplyAdd(M128D const a, M128D const b, M128D const c) {
  static_assert(fma_presence != FMAPresence::Unknown);
  if constexpr (fma_presence == FMAPresence::Present) {
    return FusedNegatedMultiplyAdd(a, b, c);
  } else {
    return c - a * b;
  }
}

// See [SZ05], section 2.1 for the validity of this function.
template<FMAPresence fma_presence>
DoublePrecision<M128D> TwoProductAdd(M128D const a,
                                     M128D const b,
                                     M128D const c) {
  static_assert(fma_presence != FMAPresence::Unknown);
  // Somehow `if constexpr` loses a cycle on MSVC 17.
  if (fma_presence == FMAPresence::Present) {
    DoublePrecision<M128D> result(uninitialized);
    result.value = FusedMultiplyAdd(a, b, c);
    result.error = FusedMultiplySubtract(a, b, result.value - c);
    return result;
  } else {
    auto result = VeltkampDekkerProduct(a, b);
    result += c;
    return result;
  }
}

// See [SZ05], section 2.1 for the validity of this function.
template<FMAPresence fma_presence>
DoublePrecision<M128D> TwoProductNegatedAdd(M128D const a,
                                            M128D const b,
                                            M128D const c) {
  static_assert(fma_presence != FMAPresence::Unknown);
  // Somehow `if constexpr` loses a cycle on MSVC 17.
  if (fma_presence == FMAPresence::Present) {
    DoublePrecision<M128D> result(uninitialized);
    result.value = FusedNegatedMultiplyAdd(a, b, c);
    result.error = FusedNegatedMultiplySubtract(a, b, result.value - c);
    return result;
  } else {
    auto result = VeltkampDekkerProduct(-a, b);
    result += c;
    return result;
  }
}

inline std::int64_t AccurateTableIndex(Argument const abs_x) {
  // This function computes the index in the accurate table:
  // 1. A suitable (large) power of 2 is added to the argument so that the last
  //    bit of the mantissa of the result corresponds to units of 1/512 of the
  //    argument.  As part of this addition, an interval of radius 1/1024 around
  //    an integral multiple of 1/512 is correctly rounded to that integral
  //    multiple.
  // 2. An `and` operation is used to only retain the last 9 bits of the
  //    mantissa.
  // 3. The result is interpreted as an integer and returned as the index.
  return (m128d::mantissa_index_bits &
          (abs_x + m128d::accurate_table_index_addend))
      .Bits<std::int64_t>();
}

// Evaluates the sum `y + δy` and performs the rounding test using the technique
// described in [Mul+10], section 11.6.3.  If rounding is safe, returns the sum;
// otherwise, returns `NaN`.  `y` is always positive.  `δy` may be positive or
// negative.
template<FMAPresence fma_presence, double e>
Value DetectDangerousRounding(Value const y, Value const δy) {
  // We don't check that `δy` is not NaN because that's how we trigger fallback
  // to the slow path.
  DCHECK(y == y);
  DoublePrecision<M128D> const sum = QuickTwoSum(y, δy);
  auto const& value = sum.value;
  auto const& error = sum.error;
  auto const muller_test_expression =
      MaybeFusedMultiplyAdd<fma_presence>(error, M128D(e), value);
  OSACA_IF(value == muller_test_expression) {
    return value;
  } else {
#if _DEBUG
    LOG_IF(ERROR, value == value && error == error)
        << std::setprecision(25) << y << " " << std::hexfloat << value << " "
        << error << " " << e;
#endif
    return m128d::quiet_NaN;
  }
}

template<FMAPresence fma_presence, bool preserve_sign>
FORCE_INLINE(inline)
void Reduce(Argument const x,
            DoublePrecision<Argument>& x_reduced,
            std::int64_t& quadrant) {
  Argument const abs_x = Abs(x);
  OSACA_IF(abs_x < π / 4) {
    x_reduced.value = x;
    x_reduced.error = m128d::zero;
    quadrant = 0;
    return;
  } OSACA_ELSE_IF(abs_x <= two_term_x_threshold) {
    // We are not very sensitive to rounding errors in this expression, because
    // in the worst case it could cause the reduced angle to jump from the
    // vicinity of π / 4 to the vicinity of -π / 4 with appropriate adjustment
    // of the quadrant.
    M128D const sign = Sign(x);
    M128D n_double =
        MaybeFusedMultiplyAdd<fma_presence>(
            abs_x, m128d::two_over_π, m128d::mantissa_reduce_shifter) -
        m128d::mantissa_reduce_shifter;

    // Don't move the computation of `n` after the if, it generates some extra
    // moves.
    Argument y;
    std::int64_t n;
    if constexpr (preserve_sign) {
      n_double = n_double ^ sign;
      n = _mm_cvtsd_si64(static_cast<__m128d>(n_double));
      y = MaybeFusedNegatedMultiplyAdd<fma_presence>(n_double, m128d::C₁, x);
    } else {
      n = _mm_cvtsd_si64(static_cast<__m128d>(n_double));
      y = MaybeFusedNegatedMultiplyAdd<fma_presence>(n_double, m128d::C₁, abs_x);
    }

    Argument const δy = n_double * m128d::δC₁;
    x_reduced = TwoDifference(y, δy);
    Argument const abs_x_reduced_value = Abs(x_reduced.value);
    OSACA_IF(abs_x_reduced_value >= two_term_x_reduced_threshold) {
      quadrant = n & 0b11;
      return;
    }
  } OSACA_ELSE_IF(abs_x <= three_term_x_threshold) {
    // Same code as above.
    M128D const sign = Sign(x);
    M128D n_double =
        MaybeFusedMultiplyAdd<fma_presence>(
            abs_x, m128d::two_over_π, m128d::mantissa_reduce_shifter) -
        m128d::mantissa_reduce_shifter;

    Argument y;
    std::int64_t n;
    if constexpr (preserve_sign) {
      n_double = n_double ^ sign;
      n = _mm_cvtsd_si64(static_cast<__m128d>(n_double));
      y = MaybeFusedNegatedMultiplyAdd<fma_presence>(n_double, m128d::C₂, x);
    } else {
      n = _mm_cvtsd_si64(static_cast<__m128d>(n_double));
      y = MaybeFusedNegatedMultiplyAdd<fma_presence>(n_double, m128d::C₂, abs_x);
    }

    Argument const yʹ = n_double * m128d::Cʹ₂;
    Argument const δy = n_double * m128d::δC₂;
    auto const z = QuickTwoSum(yʹ, δy);
    x_reduced = y - z;
    Argument const abs_x_reduced_value = Abs(x_reduced.value);
    OSACA_IF(abs_x_reduced_value >= three_term_x_reduced_threshold) {
      quadrant = n & 0b11;
      return;
    }
  }
  x_reduced.value = m128d::zero;
  x_reduced.error = m128d::quiet_NaN;
}

template<FMAPresence fma_presence>
Value SinPolynomial(Argument const x) {
  using Polynomial1 =
      HornerEvaluator<Value, Argument, 1, FMAPolicy::Auto, fma_presence>;
  return Polynomial1::Evaluate({m128d::sin_0, m128d::sin_1}, x);
}

template<FMAPresence fma_presence>
Value SinPolynomialNearZero(Argument const x) {
  using Polynomial1 =
      HornerEvaluator<Value, Argument, 1, FMAPolicy::Auto, fma_presence>;
  return Polynomial1::Evaluate(
      {m128d::sin_near_zero_0, m128d::sin_near_zero_1}, x);
}

template<FMAPresence fma_presence>
Value CosPolynomial(Argument const x) {
  using Polynomial1 =
      HornerEvaluator<Value, Argument, 1, FMAPolicy::Auto, fma_presence>;
  return Polynomial1::Evaluate({m128d::cos_0, m128d::cos_1}, x);
}

template<FMAPresence fma_presence>
FORCE_INLINE(inline)
Value SinImplementation(DoublePrecision<Argument> const x_reduced) {
  auto const x̃ = x_reduced.value;
  auto const δx̃ = x_reduced.error;
  auto const abs_x̃ = Abs(x̃);
  OSACA_IF(abs_x̃ < sin_near_zero_cutoff) {
    auto const x̃² = x̃ * x̃;
    auto const x̃³ = x̃² * x̃;
    auto const x̃³_term = MaybeFusedMultiplyAdd<fma_presence>(
        x̃³, SinPolynomialNearZero<fma_presence>(x̃²), δx̃);
    return DetectDangerousRounding<fma_presence, sin_near_zero_e>(
        x̃, x̃³_term);
  } else {
    auto const sign = Sign(x̃);
    auto const k = AccurateTableIndex(abs_x̃);
    auto const& accurate_values = SinCosAccurateTable[k];
    M128D const xₖ(accurate_values.x);
    M128D const sin_xₖ(accurate_values.sin_x);
    M128D const cos_xₖ(accurate_values.cos_x);
    // [GB91] incorporates `δx̃` in the computation of `h`.  However, `x̃` and
    // `δx̃` don't overlap and in the first interval `x̃` and `h` may be of the
    // same order of magnitude.  Instead we incorporate the terms in `δx̃` and
    // `δx̃ * h` later in the computation.  Note that the terms in `δx̃ * h²` and
    // higher are *not* computed because they don't matter.
    auto const h = abs_x̃ - xₖ;

    // The sign of the argument must be applied to the result.  It's best to do
    // this by applying it to elements of the computation that are available
    // early.
    M128D const signed_sin_xₖ = sign ^ sin_xₖ;
    M128D const signed_cos_xₖ = sign ^ cos_xₖ;
    M128D const signed_δx̃ = sign ^ δx̃;

    DoublePrecision<M128D> const sin_xₖ_plus_h_cos_xₖ =
        TwoProductAdd<fma_presence>(signed_cos_xₖ, h, signed_sin_xₖ);
    auto const h² = h * h;
    auto const h³ = h² * h;
    auto const h_plus_δx̃² = h * ((signed_δx̃ + signed_δx̃) + h);
    auto const polynomial_term =
        MaybeFusedMultiplyAdd<fma_presence>(
            signed_cos_xₖ,
            MaybeFusedMultiplyAdd<fma_presence>(
                h³, SinPolynomial<fma_presence>(h²), signed_δx̃),
            (signed_sin_xₖ * h_plus_δx̃²) * CosPolynomial<fma_presence>(h²)) +
        sin_xₖ_plus_h_cos_xₖ.error;
    return DetectDangerousRounding<fma_presence, sin_e>(
        sin_xₖ_plus_h_cos_xₖ.value, polynomial_term);
  }
}

template<FMAPresence fma_presence>
FORCE_INLINE(inline)
Value CosImplementation(DoublePrecision<Argument> const x_reduced) {
  auto const x̃ = x_reduced.value;
  auto const δx̃ = x_reduced.error;
  auto const abs_x̃ = Abs(x̃);
  auto const sign = Sign(x̃);
  auto const signed_δx̃ = sign ^ δx̃;
  auto const k = AccurateTableIndex(abs_x̃);
  auto const& accurate_values = SinCosAccurateTable[k];
  M128D const xₖ(accurate_values.x);
  M128D const sin_xₖ(accurate_values.sin_x);
  M128D const cos_xₖ(accurate_values.cos_x);
  // [GB91] incorporates `δx̃` in the computation of `h`.  However, `x̃` and `δx̃`
  // don't overlap and in the first interval `x̃` and `h` may be of the same
  // order of magnitude.  Instead we incorporate the terms in `δx̃` and `δx̃ * h`
  // later in the computation.  Note that the terms in `δx̃ * h²` and higher are
  // *not* computed because they don't matter.
  auto const h = abs_x̃ - xₖ;

  DoublePrecision<M128D> const cos_xₖ_minus_h_sin_xₖ =
      TwoProductNegatedAdd<fma_presence>(sin_xₖ, h, cos_xₖ);
  auto const h² = h * h;
  auto const h³ = h² * h;
  auto const h_plus_δx̃² = h * ((signed_δx̃ + signed_δx̃) + h);
  auto const polynomial_term =
      MaybeFusedNegatedMultiplyAdd<fma_presence>(
          sin_xₖ,
          MaybeFusedMultiplyAdd<fma_presence>(
              h³, SinPolynomial<fma_presence>(h²), signed_δx̃),
          (cos_xₖ * h_plus_δx̃²) * CosPolynomial<fma_presence>(h²)) +
      cos_xₖ_minus_h_sin_xₖ.error;
  return DetectDangerousRounding<fma_presence, cos_e>(
      cos_xₖ_minus_h_sin_xₖ.value, polynomial_term);
}

template<FMAPresence fma_presence>
FORCE_INLINE(inline)
SC<Value> SinCosImplementation(DoublePrecision<Argument> const x_reduced) {
  SC<Value> m128ds;
  auto const x̃ = x_reduced.value;
  auto const δx̃ = x_reduced.error;
  auto const abs_x̃ = Abs(x̃);
  auto const sign = Sign(x̃);
  auto const k = AccurateTableIndex(abs_x̃);
  auto const& accurate_values = SinCosAccurateTable[k];
  M128D const xₖ(accurate_values.x);
  M128D const sin_xₖ(accurate_values.sin_x);
  M128D const cos_xₖ(accurate_values.cos_x);
  // [GB91] incorporates `δx̃` in the computation of `h`.  However, `x̃` and `δx̃`
  // don't overlap and in the first interval `x̃` and `h` may be of the same
  // order of magnitude.  Instead we incorporate the terms in `δx̃` and `δx̃ * h`
  // later in the computation.  Note that the terms in `δx̃ * h²` and higher are
  // *not* computed because they don't matter.
  auto const h = abs_x̃ - xₖ;

  // The sign of the argument must be applied to the result.  It's best to do
  // this by applying it to elements of the computation that are available
  // early.
  M128D const signed_sin_xₖ = sign ^ sin_xₖ;
  M128D const signed_cos_xₖ = sign ^ cos_xₖ;
  M128D const signed_δx̃ = sign ^ δx̃;

  DoublePrecision<M128D> const sin_xₖ_plus_h_cos_xₖ =
      TwoProductAdd<fma_presence>(signed_cos_xₖ, h, signed_sin_xₖ);
  DoublePrecision<M128D> const cos_xₖ_minus_h_sin_xₖ =
      TwoProductNegatedAdd<fma_presence>(sin_xₖ, h, cos_xₖ);
  auto const h² = h * h;
  auto const h³ = h² * h;
  auto const h_plus_δx̃² = h * ((signed_δx̃ + signed_δx̃) + h);

  auto const h³_sin_polynomial = MaybeFusedMultiplyAdd<fma_presence>(
      h³, SinPolynomial<fma_presence>(h²), signed_δx̃);
  auto const h_plus_e²_cos_polynomial =
      h_plus_δx̃² * CosPolynomial<fma_presence>(h²);

  auto const sin_polynomial_term =
      MaybeFusedMultiplyAdd<fma_presence>(
          signed_cos_xₖ,
          h³_sin_polynomial,
          signed_sin_xₖ * h_plus_e²_cos_polynomial) +
      sin_xₖ_plus_h_cos_xₖ.error;
  auto const cos_polynomial_term =
      MaybeFusedNegatedMultiplyAdd<fma_presence>(
          sin_xₖ,
          h³_sin_polynomial,
          cos_xₖ * h_plus_e²_cos_polynomial) +
      cos_xₖ_minus_h_sin_xₖ.error;
  m128ds.cos = DetectDangerousRounding<fma_presence, cos_e>(
      cos_xₖ_minus_h_sin_xₖ.value, cos_polynomial_term);
  OSACA_IF(abs_x̃ < sin_near_zero_cutoff) {
    auto const x̃² = x̃ * x̃;
    auto const x̃³ = x̃² * x̃;
    auto const x̃³_term = MaybeFusedMultiplyAdd<fma_presence>(
        x̃³, SinPolynomialNearZero<fma_presence>(x̃²), δx̃);
    m128ds.sin =
        DetectDangerousRounding<fma_presence, sin_near_zero_e>(x̃, x̃³_term);
  } else {
    m128ds.sin = DetectDangerousRounding<fma_presence, sin_e>(
        sin_xₖ_plus_h_cos_xₖ.value, sin_polynomial_term);
  }
  return m128ds;
}

template<FMAPresence fma_presence>
double __cdecl Sin(double x) {
  OSACA_FUNCTION_BEGIN(x, <fma_presence>);
  DoublePrecision<Argument> x_reduced;
  std::int64_t quadrant;
  double value;
  Reduce<fma_presence, /*preserve_sign=*/true>(
      M128D(x), x_reduced, quadrant);
  OSACA_IF(quadrant & 0b1) {
    value = static_cast<double>(CosImplementation<fma_presence>(x_reduced));
  } else {
    value = static_cast<double>(SinImplementation<fma_presence>(x_reduced));
  }
  OSACA_IF(value != value) {
    if (slow_path_sin_callback != nullptr) {
      slow_path_sin_callback(x);
    }
    OSACA_RETURN(cr_sin(x));
  } OSACA_ELSE_IF(quadrant & 0b10) {
    OSACA_RETURN(-value);
  } else {
    OSACA_RETURN(value);
  }
}

template<FMAPresence fma_presence>
double __cdecl Cos(double x) {
  OSACA_FUNCTION_BEGIN(x, <fma_presence>);
  DoublePrecision<Argument> x_reduced;
  std::int64_t quadrant;
  double value;
  Reduce<fma_presence, /*preserve_sign=*/false>(
      M128D(x), x_reduced, quadrant);
  OSACA_IF(quadrant & 0b1) {
    value = static_cast<double>(SinImplementation<fma_presence>(x_reduced));
  } else {
    value = static_cast<double>(CosImplementation<fma_presence>(x_reduced));
  }
  OSACA_IF(value != value) {
    if (slow_path_cos_callback != nullptr) {
      slow_path_cos_callback(x);
    }
    OSACA_RETURN(cr_cos(x));
  } OSACA_ELSE_IF(quadrant == 1 || quadrant == 2) {
    OSACA_RETURN(-value);
  } else {
    OSACA_RETURN(value);
  }
}

template<FMAPresence fma_presence>
SC<double> __cdecl SinCos(double x) {
  OSACA_FUNCTION_BEGIN(x, <fma_presence>);
  SC<double> values;
  DoublePrecision<Argument> x_reduced;
  std::int64_t quadrant;
  Reduce<fma_presence, /*preserve_sign=*/true>(
      M128D(x), x_reduced, quadrant);
  auto m128ds = SinCosImplementation<fma_presence>(x_reduced);
  OSACA_IF(quadrant & 0b1) {
    values.sin = static_cast<double>(m128ds.cos);
    values.cos = static_cast<double>(m128ds.sin);
  } else {
    values.sin = static_cast<double>(m128ds.sin);
    values.cos = static_cast<double>(m128ds.cos);
  }
  OSACA_IF(values.sin != values.sin || values.cos != values.cos) {
    if (slow_path_sin_callback != nullptr) {
      slow_path_sin_callback(x);
    }
    if (slow_path_cos_callback != nullptr) {
      slow_path_cos_callback(x);
    }
    OSACA_RETURN((SC<double>{.sin = cr_sin(x), .cos = cr_cos(x)}));
  }
  OSACA_IF(quadrant & 0b10) {
    values.sin = -values.sin;
  }
  OSACA_IF(quadrant == 1 || quadrant == 2) {
    values.cos = -values.cos;
  }
  OSACA_RETURN(values);
}

void SetSlowPathsCallbacks(SlowPathCallback sin_cb, SlowPathCallback cos_cb) {
  slow_path_sin_callback = std::move(sin_cb);
  slow_path_cos_callback = std::move(cos_cb);
}

template double __cdecl Sin<FMAPresence::Absent>(double x);
template double __cdecl Sin<FMAPresence::Present>(double x);
template double __cdecl Cos<FMAPresence::Absent>(double x);
template double __cdecl Cos<FMAPresence::Present>(double x);
template SC<double> __cdecl SinCos<FMAPresence::Absent>(double x);
template SC<double> __cdecl SinCos<FMAPresence::Present>(double x);

}  // namespace internal
}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
