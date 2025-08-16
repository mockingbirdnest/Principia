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
// TODO(phl): Update the code to match the notation of the document.
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
    /* Not NaN is the only part that matters; used at the end of the    */   \
    /* top-level functions to determine whether to call the slow path.  */   \
    constexpr double value = 1;                                              \
    constexpr double muller_test_expression = value;                         \
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
constexpr double two_term_θ_threshold = π / 2 * (1LL << κ₁);
constexpr double three_term_θ_threshold = π / 2 * (1LL << κ₂);
constexpr double C₁ = 0x1.921F'B544'42D0'0p0;
constexpr double δC₁ = 0x1.8469'898C'C517'0p-48;
constexpr double C₂ = 0x1.921F'B544'4000'0p0;
constexpr double Cʹ₂ = 0x1.68C2'34C4'C000'0p-39;
constexpr double δC₂ = 0x1.98A2'E037'0734'5p-77;
constexpr double two_term_θ_reduced_threshold =
    1.0 / (1LL << (-(κ₁ + κʹ₁ + κ₃ - std::numeric_limits<double>::digits + 2)));
constexpr double three_term_θ_reduced_threshold =
    (1.0 / (1LL << (-(κ₃ - std::numeric_limits<double>::digits)))) *
    ((1LL << (-(κ₂ + κʹ₂ + κʺ₂ - std::numeric_limits<double>::digits + 2))) +
     4);

M128D const mantissa_reduce_shifter(
    static_cast<double>(1LL << (std::numeric_limits<double>::digits - 1)));
M128D const accurate_table_index_addend(static_cast<double>(
    1LL << (std::numeric_limits<double>::digits - table_spacing_bits - 1)));

constexpr double sin_near_zero_e = 0x1.0000'AD82'A723'6p0;  // 2^-70.561.
constexpr double sin_e = 0x1.0002'6013'6BD9'Dp0;  // 2^-68.751.
constexpr double cos_e = 0x1.0001'B836'988A'Dp0;  // 2^-69.217.

template<FMAPolicy fma_policy>
using Polynomial1 = HornerEvaluator<Value, Argument, 1, fma_policy>;

// Pointers used for indirect calls, set by `StaticInitialization`.
double (__cdecl *cos)(double θ) = nullptr;
double (__cdecl *sin)(double θ) = nullptr;

// Forward declarations needed by the OSACA macros.
template<FMAPolicy fma_policy>
double __cdecl Sin(double θ);
template<FMAPolicy fma_policy>
double __cdecl Cos(double θ);

namespace masks {
M128D const mantissa_index_bits = M128D::MakeFromBits(0x0000'0000'0000'01ffull);
}  // namespace masks

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
  return (masks::mantissa_index_bits & (abs_x + accurate_table_index_addend))
      .Bits<std::int64_t>();
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

M128D const nan(std::numeric_limits<double>::quiet_NaN());

// Evaluates the sum `x + Δx` and performs the rounding test using the technique
// described in [Mul+10], section 11.6.3.  If rounding is safe, returns the sum;
// otherwise, returns `NaN`.  `x` is always positive.  `Δx` may be positive or
// negative.
template<FMAPolicy fma_policy, double e>
Value DetectDangerousRounding(Value const x, Value const Δx) {
  // We don't check that `Δx` is not NaN because that's how we trigger fallback
  // to the slow path.
  DCHECK(x == x);
  DoublePrecision<M128D> const sum = QuickTwoSum(x, Δx);
  auto const& value = sum.value;
  auto const& error = sum.error;
  auto const muller_test_expression =
      FusedMultiplyAdd<fma_policy>(error, M128D(e), value);
  OSACA_IF(value == muller_test_expression) {
    return value;
  } else {
#if _DEBUG
    LOG_IF(ERROR, value == value && error == error)
        << std::setprecision(25) << x << " " << std::hexfloat << value << " "
        << error << " " << e;
#endif
    return nan;
  }
}

M128D const zero(0.0);
M128D const two_over_π(2.0 / π);
M128D const C1(C₁);
M128D const δC1(δC₁);
M128D const C2(C₂);
M128D const Cp2(Cʹ₂);
M128D const δC2(δC₂);

template<FMAPolicy fma_policy, bool preserve_sign>
FORCE_INLINE(inline)
void Reduce(Argument const θ,
            DoublePrecision<Argument>& θ_reduced,
            std::int64_t& quadrant) {
  Argument const abs_θ = Abs(θ);
  OSACA_IF(abs_θ < π / 4) {
    θ_reduced.value = θ;
    θ_reduced.error = zero;
    quadrant = 0;
    return;
  } OSACA_ELSE_IF(abs_θ <= two_term_θ_threshold) {
    // We are not very sensitive to rounding errors in this expression, because
    // in the worst case it could cause the reduced angle to jump from the
    // vicinity of π / 4 to the vicinity of -π / 4 with appropriate adjustment
    // of the quadrant.
    M128D const sign = Sign(θ);
    M128D n_double = FusedMultiplyAdd<fma_policy>(
                         abs_θ, two_over_π, mantissa_reduce_shifter) -
                     mantissa_reduce_shifter;

    // Don't move the computation of `n` after the if, it generates some extra
    // moves.
    Argument y;
    std::int64_t n;
    if constexpr (preserve_sign) {
      n_double = n_double ^ sign;
      n = _mm_cvtsd_si64(static_cast<__m128d>(n_double));
      y = FusedNegatedMultiplyAdd<fma_policy>(n_double, C1, θ);
    } else {
      n = _mm_cvtsd_si64(static_cast<__m128d>(n_double));
      y = FusedNegatedMultiplyAdd<fma_policy>(n_double, C1, abs_θ);
    }

    Argument const δy = n_double * δC1;
    θ_reduced = TwoDifference(y, δy);
    OSACA_IF(θ_reduced.value <= -two_term_θ_reduced_threshold ||
             θ_reduced.value >= two_term_θ_reduced_threshold) {
      quadrant = n & 0b11;
      return;
    }
  } OSACA_ELSE_IF(abs_θ <= three_term_θ_threshold) {
    // Same code as above.
    M128D const sign = Sign(θ);
    M128D n_double = FusedMultiplyAdd<fma_policy>(
                         abs_θ, two_over_π, mantissa_reduce_shifter) -
                     mantissa_reduce_shifter;

    Argument y;
    std::int64_t n;
    if constexpr (preserve_sign) {
      n_double = n_double ^ sign;
      n = _mm_cvtsd_si64(static_cast<__m128d>(n_double));
      y = FusedNegatedMultiplyAdd<fma_policy>(n_double, C2, θ);
    } else {
      n = _mm_cvtsd_si64(static_cast<__m128d>(n_double));
      y = FusedNegatedMultiplyAdd<fma_policy>(n_double, C2, abs_θ);
    }

    Argument const yʹ = n_double * Cp2;
    Argument const δy = n_double * δC2;
    auto const z = QuickTwoSum(yʹ, δy);
    θ_reduced = y - z;
    OSACA_IF(θ_reduced.value <= -three_term_θ_reduced_threshold ||
             θ_reduced.value >= three_term_θ_reduced_threshold) {
      quadrant = n & 0b11;
      return;
    }
  }
  θ_reduced.value = zero;
  θ_reduced.error = nan;
}

M128D const s0(-0x1.5555'5555'5555'5p-3);
M128D const s1(0x1.1111'10A8'20AE'Cp-7);

template<FMAPolicy fma_policy>
Value SinPolynomial(Argument const x) {
  return Polynomial1<fma_policy>::Evaluate({s0, s1}, x);
}

M128D const s00(-0x1.5555'5555'5555'5p-3);
M128D const s01(0x1.1111'10B4'0E88'Ap-7);

template<FMAPolicy fma_policy>
Value SinPolynomialNearZero(Argument const x) {
  return Polynomial1<fma_policy>::Evaluate({s00, s01}, x);
}

M128D const c0(-0x1.FFFF'FFFF'FFFF'Dp-2);
M128D const c1(0x1.5555'549D'B0A9'5p-5);

template<FMAPolicy fma_policy>
Value CosPolynomial(Argument const x) {
  return Polynomial1<fma_policy>::Evaluate({c0, c1}, x);
}

M128D const two(2.0);

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value SinImplementation(DoublePrecision<Argument> const θ_reduced) {
  auto const x = θ_reduced.value;
  auto const e = θ_reduced.error;
  auto const abs_x = Abs(x);
  OSACA_IF(abs_x < sin_near_zero_cutoff) {
    auto const x² = x * x;
    auto const x³ = x² * x;
    auto const x³_term = FusedMultiplyAdd<fma_policy>(
        x³, SinPolynomialNearZero<fma_policy>(x²), e);
    return DetectDangerousRounding<fma_policy, sin_near_zero_e>(x, x³_term);
  } else {
    auto const sign = Sign(x);
    auto const e_abs = e ^ sign;
    auto const i = AccurateTableIndex(abs_x);
    auto const& accurate_values = SinCosAccurateTable[i];
    M128D const x₀(accurate_values.x);
    M128D const sin_x₀(accurate_values.sin_x);
    M128D const cos_x₀(accurate_values.cos_x);
    // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
    // don't overlap and in the first interval `x` and `h` may be of the same
    // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
    // later in the computation.  Note that the terms in `e * h²` and higher are
    // *not* computed because they don't matter.
    auto const h = abs_x - x₀;

    DoublePrecision<M128D> const sin_x₀_plus_h_cos_x₀ =
        TwoProductAdd<fma_policy>(cos_x₀, h, sin_x₀);
    auto const h² = h * h;
    auto const h³ = h² * h;
    auto const h_plus_e_abs² =
        h * FusedMultiplyAdd<fma_policy>(two, e_abs, h);
    auto const polynomial_term =
        FusedMultiplyAdd<fma_policy>(
            cos_x₀,
            FusedMultiplyAdd<fma_policy>(
                h³, SinPolynomial<fma_policy>(h²), e_abs),
            (sin_x₀ * h_plus_e_abs²) * CosPolynomial<fma_policy>(h²)) +
        sin_x₀_plus_h_cos_x₀.error;
    return sign ^ DetectDangerousRounding<fma_policy, sin_e>(
                      sin_x₀_plus_h_cos_x₀.value, polynomial_term);
  }
}

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value CosImplementation(DoublePrecision<Argument> const θ_reduced) {
  auto const x = θ_reduced.value;
  auto const e = θ_reduced.error;
  auto const abs_x = Abs(x);
  auto const sign = Sign(x);
  auto const e_abs = e ^ sign;
  auto const i = AccurateTableIndex(abs_x);
  auto const& accurate_values = SinCosAccurateTable[i];
  M128D const x₀(accurate_values.x);
  M128D const sin_x₀(accurate_values.sin_x);
  M128D const cos_x₀(accurate_values.cos_x);
  // [GB91] incorporates `e` in the computation of `h`.  However, `x` and `e`
  // don't overlap and in the first interval `x` and `h` may be of the same
  // order of magnitude.  Instead we incorporate the terms in `e` and `e * h`
  // later in the computation.  Note that the terms in `e * h²` and higher are
  // *not* computed because they don't matter.
  auto const h = abs_x - x₀;

  DoublePrecision<M128D> const cos_x₀_minus_h_sin_x₀ =
      TwoProductNegatedAdd<fma_policy>(sin_x₀, h, cos_x₀);
  auto const h² = h * h;
  auto const h³ = h² * h;
  auto const h_plus_e_abs² =
      h * FusedMultiplyAdd<fma_policy>(two, e_abs, h);
  // TODO(phl): Redo the error analysis.
  auto const polynomial_term =
      FusedNegatedMultiplyAdd<fma_policy>(
          sin_x₀,
          FusedMultiplyAdd<fma_policy>(
              h³, SinPolynomial<fma_policy>(h²), e_abs),
          h_plus_e_abs² * (cos_x₀ * CosPolynomial<fma_policy>(h²))) +
      cos_x₀_minus_h_sin_x₀.error;
  return DetectDangerousRounding<fma_policy, cos_e>(cos_x₀_minus_h_sin_x₀.value,
                                                    polynomial_term);
}

template<FMAPolicy fma_policy>
double __cdecl Sin(double θ) {
  OSACA_FUNCTION_BEGIN(θ, <fma_policy>);
  DoublePrecision<Argument> θ_reduced;
  std::int64_t quadrant;
  double value;
  Reduce<fma_policy, /*preserve_sign=*/true>(M128D(θ), θ_reduced, quadrant);
  OSACA_IF(quadrant & 0b1) {
    value = static_cast<double>(CosImplementation<fma_policy>(θ_reduced));
  } else {
    value = static_cast<double>(SinImplementation<fma_policy>(θ_reduced));
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
double __cdecl Cos(double θ) {
  OSACA_FUNCTION_BEGIN(θ, <fma_policy>);
  DoublePrecision<Argument> θ_reduced;
  std::int64_t quadrant;
  double value;
  Reduce<fma_policy, /*preserve_sign=*/false>(M128D(θ), θ_reduced, quadrant);
  OSACA_IF(quadrant & 0b1) {
    value = static_cast<double>(SinImplementation<fma_policy>(θ_reduced));
  } else {
    value = static_cast<double>(CosImplementation<fma_policy>(θ_reduced));
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

double __cdecl Sin(double const θ) {
  return sin(θ);
}

double __cdecl Cos(double const θ) {
  return cos(θ);
}

}  // namespace internal
}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
