#include "numerics/cbrt.hpp"

#include <pmmintrin.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <utility>

#include "glog/logging.h"
#include "numerics/double_precision.hpp"
#include "numerics/fma.hpp"
#include "quantities/elementary_functions.hpp"


#define OSACA_ANALYSED_FUNCTION Cbrt
#define UNDER_OSACA_HYPOTHESES(expression)                                 \
  [&] {                                                                     \
    constexpr double y = 3;                                                \
    constexpr double abs_y = y == 0 ? 0 : y > 0 ? y : -y;                  \
    constexpr double r₀ = std::numeric_limits<double>::quiet_NaN();        \
    constexpr double r₁ = std::numeric_limits<double>::quiet_NaN();        \
    constexpr double r̃ = std::numeric_limits<double>::quiet_NaN();         \
    constexpr auto CorrectionPossiblyNeeded = [](double const r₀,          \
                                                 double const r₁,          \
                                                 double const r̃,           \
                                                 double const τ) -> bool { \
      return false;                                                        \
    };                                                                     \
    return expression;                                                     \
  }()

#if defined(OSACA_ANALYSED_FUNCTION)
#define PRINCIPIA_USE_OSACA !PRINCIPIA_MACRO_IS_EMPTY(OSACA_ANALYSED_FUNCTION)
#endif

#if PRINCIPIA_USE_OSACA

#include "intel/iacaMarks.h"

// The macros OSACA_FUNCTION_BEGIN and OSACA_RETURN are used to analyse the
// latency of a double -> double function, as measured, e.g., by the
// nanobenchmarks; this notionally corresponds to the duration of an iteration
// of a loop `x = f(x)`.
// The latency-critical path of the function is reported as the loop-carried
// dependency by OSACA, and as the critical path by IACA in throughput analysis
// mode.
// OSACA and IACA are only suitable for benchmarking a single path.  Any if
// statements, including in inline functions called by the function being
// analysed, should be replaced with OSACA_IF.  The path should be determined by
// providing any necessary constexpr expressions in UNDER_OSACA_HYPOTHESES.
// If OSACA_EVALUATE_CONDITIONS is set to 1, code will be generated to evaluate
// the conditions for the branches taken (outside the loop-carried path, as they
// would be with correct branch prediction).
// OSACA_EVALUATE_CONDITIONS can be set to 0 to get a more streamlined
// dependency graph when the conditions are irrelevant.  Note that this can
// result in artificially low port pressures.
#if 0
// Usage:
  double f(double const x) {
    double const abs_x = std::abs(x);
    if (abs_x < 0.1) {
      return x;
    } else if (x < 0) {
      return -f_positive(abs_x);
    } else {
      return f_positive(abs_x);
    }
  }
  double f_positive(double const abs_x) {
    if (abs_x > 3) {
      return 1;
    } else {
      // …
    }
  }
// becomes:
  double f(double x) {  // The argument cannot be const.
    OSACA_FUNCTION_BEGIN(x);
    double const abs_x = std::abs(x);
    OSACA_IF(abs_x < 0.1) {
      OSACA_RETURN(x);
    } OSACA_ELSE_IF(x < 0) {  // `else OSACA_IF` works, but upsets the linter.
      OSACA_RETURN(-f_positive(abs_x));
    } else {
      OSACA_RETURN(f_positive(abs_x));
    }
  }
  // Other functions can have const arguments and return normally, but need to
  // use OSACA_IF:
  double f_positive(double const abs_x) {
    OSACA_IF(abs_x > 3) {
      return 1;
    } else {
      // …
    }
  }
// To analyse it near x = 5:
#define OSACA_ANALYSED_FUNCTION f
#define UNDER_OSACA_HYPOTHESES(expression)                                 \
  [&] {                                                                    \
    constexpr double x = 5;                                                \
    /* To avoid inconsistent definitions, constexpr code can be used as */ \
    /* needed.                                                          */ \
    constexpr double abs_x = x > 0 ? x : -x;                               \
    return (expression);                                                   \
  }()
#endif
// In principle, the loop-carried dependency for function latency analysis is
// achieved by wrapping the body of the function in an infinite loop, assigning
// the result to the argument at each iteration, and adding IACA markers outside
// the loop.  There are some subtleties:
// — We need to trick the compiler into believing the loop is finite, so that it
//   doesn’t optimize away the end marker or even the function.  This is
//   achieved by exiting based on the value of `OSACA_loop_terminator`.
// — Return statements may be in if statements, and there may be several of
//   them, so they cannot be the end of a loop started unconditionally.  Instead
//   we loop with goto.
// — We need to prevent the compiler from moving the start and end markers into
//   the middle of register saving and restoring code, which would mess up the
//   dependency analysis.  This is done with additional conditional gotos.
// — Some volatile reads and writes are used to clarify identity of the
//   registers in the generated code (where the names of `OSACA_result` and, if
//   `OSACA_CARRY_LOOP_THROUGH_REGISTER` is set to 0, `OSACA_loop_carry` appear
//   in movsd instructions).
//
// Carrying the loop dependency through a memory load and store can make the
// dependency graph easier to understand, as it forces any usage of the input to
// depend on the initial movsd, with the loop carried by a single backward edge
// to that initial movsd.
// If the loop is carried through a register, multiple usages of the input may
// result in multiple back edges from the final instruction that computed the
// result.  Carrying the loop through the memory could also potentially prevent
// the compiler from reusing intermediate values in the next iteration, e.g., if
// the the computation of f(x) depends on -x and produces -f(x) before f(x), as
// in an even function defined in terms of its positive half, the compiler might
// reuse -f(x₀)=-x₁ instead of computing -x₁ from x₁=f(x₀).  However:
// — it adds a spurious move to the latency;
// — some tools (IACA) cannot see the dependency through memory.
// Set OSACA_CARRY_LOOP_THROUGH_REGISTER to 0 to carry the loop through memory.

#define OSACA_EVALUATE_CONDITIONS 1
#define OSACA_CARRY_LOOP_THROUGH_REGISTER 1

static bool volatile OSACA_loop_terminator = false;

#define OSACA_FUNCTION_BEGIN(arg)                               \
  double OSACA_LOOP_CARRY_QUALIFIER OSACA_loop_carry = arg;     \
  OSACA_outer_loop:                                             \
  if constexpr (std::string_view(__func__) ==                   \
                STRINGIFY_EXPANSION(OSACA_ANALYSED_FUNCTION)) { \
    IACA_VC64_START;                                            \
  }                                                             \
  _Pragma("warning(push)");                                     \
  _Pragma("warning(disable : 4102)");                           \
  OSACA_loop:                                                   \
  _Pragma("warning(pop)");                                      \
  arg = OSACA_loop_carry

#define OSACA_RETURN(result)                                                \
  do {                                                                      \
    if constexpr (std::string_view(__func__) ==                             \
                  STRINGIFY_EXPANSION(OSACA_ANALYSED_FUNCTION)) {           \
      OSACA_loop_carry = (result);                                          \
      if (!OSACA_loop_terminator) {                                         \
        goto OSACA_loop;                                                    \
      }                                                                     \
      double volatile OSACA_result = OSACA_loop_carry;                      \
      IACA_VC64_END;                                                        \
      /* The outer loop prevents the the start and end marker from being */ \
      /* interleaved with register saving and restoring moves.           */ \
      if (!OSACA_loop_terminator) {                                         \
        goto OSACA_outer_loop;                                              \
      }                                                                     \
      return OSACA_result;                                                  \
    } else {                                                                \
      return (result);                                                      \
    }                                                                       \
  } while (false)

#if OSACA_CARRY_LOOP_THROUGH_REGISTER
#define OSACA_LOOP_CARRY_QUALIFIER
#else
#define OSACA_LOOP_CARRY_QUALIFIER volatile
#endif

// The branch not taken, determined by evaluating the condition
// `UNDER_OSACA_HYPOTHESES`, is eliminated by `if constexpr`; the condition is
// also compiled normally and assigned to a boolean; whether this results in any
// generated code depends on `OSACA_EVALUATE_CONDITIONS`.  Note that, with
// `OSACA_EVALUATE_CONDITIONS`, in  `OSACA_IF(p) { } OSACA_ELSE_IF(q) { }`, if
// `p` holds `UNDER_OSACA_HYPOTHESES`, code is generated to evaluate `p`, but
// not `q`.

#define OSACA_IF(condition)                                               \
  if constexpr (bool OSACA_CONDITION_QUALIFIER OSACA_computed_condition = \
                    (condition);                                          \
                UNDER_OSACA_HYPOTHESES(condition))

#if OSACA_EVALUATE_CONDITIONS
#define OSACA_CONDITION_QUALIFIER volatile
#else
#define OSACA_CONDITION_QUALIFIER
#endif

#else  // if !PRINCIPIA_USE_OSACA

#define OSACA_FUNCTION_BEGIN(arg)
#define OSACA_RETURN(result) return (result)
#define OSACA_IF(condition) if (condition)

#endif  // PRINCIPIA_USE_OSACA

#define OSACA_ELSE_IF else OSACA_IF  // NOLINT


namespace principia {
namespace numerics {
namespace _cbrt {
namespace internal {

using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_fma;
using namespace principia::quantities::_elementary_functions;

// The computations in this file are described in documentation/cbrt.pdf; the
// identifiers match the notation in that document.

namespace {

// See [Nie04], algorithm 10.
std::array<double, 4> NievergeltQuadruplyCompensatedStep(
    DoublePrecision<double> b,
    DoublePrecision<double> d) {
  auto const& [b₁, b₂] = b;
  auto const& [d₁, d₂] = d;
  if (std::abs(d₂) >= std::abs(b₁)) {
    std::swap(b, d);
  }
  if (std::abs(b₂) >= std::abs(d₁)) {
    double const g₄ = d₂;
    auto const [w, g₃] = TwoSum(b₂, d₁);
    auto const [g₁, g₂] = TwoSum(b₁, w);
    return {g₁, g₂, g₃, g₄};
  } else {
    auto const [w, g₄] = TwoSum(b₂, d₂);
    auto const [u, v] = TwoSum(d₁, b₁);
    auto const [x, g₃] = TwoSum(v, w);
    auto const [g₁, g₂] = TwoSum(u, x);
    return {g₁, g₂, g₃, g₄};
  }
}

// See [Nie04], algorithm 8; note that we use 0-based indices where Nievergelt
// uses 1-based indices.  In order to avoid allocations, we do not return a
// vector, but instead an array of the maximum possible size.  Should it be
// necessary to return the number of nonzero elements of the result, some sort
// of bounded vector would probably be useful.  Note that this size is k in
// Nievergelt’s procedure, but with our indexing it is k + 1 at the return
// statement.
template<std::size_t n>
std::array<double, n> PriestNievergeltNormalize(std::array<double, n> const f) {
  std::array<double, n> s{};
  int k = 0;
  s[0] = f[0];
  for (int j = 1; j < n; ++j) {
    auto const [c, d] = TwoSum(s[k], f[j]);
    s[k] = c;
    if (d != 0) {
      int l = k - 1;
      k = k + 1;
      while (l >= 0) {
        auto const [cʹ, dʹ] = TwoSum(s[l], s[l + 1]);
        s[l] = cʹ;
        if (dʹ == 0) {
          k = k - 1;
        } else {
          s[l + 1] = dʹ;
        }
        l = l - 1;
      }
      s[k] = d;
    }
  }
  return s;
}

bool CorrectionPossiblyNeeded(double const r₀,
                              double const r₁,
                              double const r̃,
                              double const τ) {
  // The order of the conjunction matters for branch predictability here.  The
  // left-hand side is almost always false: it is true only when r = r₀ + r₁ is
  // very close to a tie (for r̃ ≠ r₀) or to a representable number (for r̃ = r₀).
  // On the other hand, it is hard to predict whether r̃ = r₀, i.e., whether r is
  // closer to a representable number than to a tie.
  return std::abs(0.5 * (r̃ - r₀) - r₁) <= τ * r₀ && r̃ != r₀;
}

double CorrectLastBit(double const y, double const r₀, double const r̃) {
  double const a = std::min(r₀, r̃);
  double const b = 0.5 * std::abs(r₀ - r̃);
  return CbrtOneBit(y, a, b) ? std::max(r₀, r̃) : a;
}

}  // namespace


namespace masks {
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d round_toward_zero_17_bits =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
static const __m128d round_toward_zero_26_bits =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFFF'F800'0000));
}  // namespace masks

namespace method_3²ᴄZ5¹ {

// No overflow or underflow occurs in intermediate computations for
// y ∈ [y₁, y₂].  The limiting expression is the denominator of step 4, wherein
// the 8th power of x, which is approximately the 8th power of ∛y, appears. The
// the 9th power of ∛y, i.e., cube of y, does not overflow for y ∈ [y₁, y₂].
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double y₁ = 0x1p-340;
constexpr double σ₁ = 0x1p-227;
constexpr double σ₁⁻³ = 1 / (σ₁ * σ₁ * σ₁);
constexpr double y₂ = 0x1p341;
constexpr double σ₂ = 0x1p227;
constexpr double σ₂⁻³ = 1 / (σ₂ * σ₂ * σ₂);
static_assert(σ₁⁻³ * y₁ == y₂);
static_assert(σ₂⁻³ * y₂ == y₁);

template<Rounding rounding>
double Cbrt(double const y) {
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(masks::sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(masks::sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);

  if (y != y) {
    // The usual logic will produce a qNaN when given a NaN, but will not
    // preserve the payload and will signal overflows (q will be a nonsensical
    // large value, and q³ will overflow).  Further, the rescaling comparisons
    // will signal the invalid operation exception for quiet NaNs (although that
    // would be easy to work around using the unordered compare intrinsics).
    return y + y;
  }

  if (abs_y < y₁) {
    if (abs_y == 0) {
      return y;
    }
    return Cbrt<rounding>(y * σ₁⁻³) * σ₁;
  } else if (abs_y > y₂) {
    if (abs_y == std::numeric_limits<double>::infinity()) {
      return y;
    }
    return Cbrt<rounding>(y * σ₂⁻³) * σ₂;
  }

  // Step 1.  The constant Γʟ²ᴄ is the result of Canon optimization for step 2.
  constexpr double Γʟ²ᴄ = 0x0.199E'9760'9F63'9F90'626F'8B97'2B3A'6249'2p0;
  constexpr std::uint64_t C = 0x2A9F'775C'D8A7'5897;
  // The fixed-point number C = ⟦(2×1023 − Γʟ²ᴄ) / 3⟧ is not representable as a
  // double; frac C + 1 = ⟦2 − Γʟ²ᴄ / 3⟧ is, and has the same last place.
  constexpr double frac_C_plus_1 = (C & 0x000F'FFFF'FFFF'FFFF) * 0x1p-52 + 1;
  // By sheer luck it happens that ⟦2 − Γʟ²ᴄ / 3⟧ = ⟦2 − ⟦⟦Γʟ²ᴄ⟧ / 3⟧⟧.
  static_assert(frac_C_plus_1 == 2 - Γʟ²ᴄ / 3);
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = _mm_cvtsd_f64(_mm_castsi128_pd(_mm_cvtsi64_si128(Q)));

  // Step 2, Lagny’s irrational method with Canon optimization, with the
  // evaluation strategy described in appendix D.
  double const q² = q * q;
  double const q⁴ = q² * q²;
  double const ξ =
      (0x1.BBA02BAFEA9B7p0 * q² + Sqrt(0x1.0030F1F8A11DAp2 * abs_y * q - q⁴)) *
      (0x1.2774CDF81A35p-2 / q);

  // Step 3.
  double const x = _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_sd(ξ), masks::round_toward_zero_17_bits));

  // Step 4, the Lagny–Schröder rational method of order 5.
  double const x³ = x * x * x;  // Exact.
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const x²_sign_y = x_sign_y * x;
  double const numerator = (x³ - abs_y) * ((10 * x³ + 16 * abs_y) * x³ + y²);
  double const denominator =
      x²_sign_y * ((15 * x³ + 51 * abs_y) * x³ + 15 * y²);
  double const Δ = numerator / denominator;
  double const r₀ = x_sign_y - Δ;
  double const r₁ = x_sign_y - r₀ - Δ;

  double const r̃ = r₀ + 2 * r₁;
  if (rounding == Rounding::Correct &&
      CorrectionPossiblyNeeded(r₀, r₁, r̃, /*τ=*/0x1.7C73DBBD9FA60p-66)) {
    return _mm_cvtsd_f64(_mm_or_pd(
        _mm_set_sd(CorrectLastBit(abs_y, std::abs(r₀), std::abs(r̃))), sign));
  }
  return r₀;
}
template double Cbrt<Rounding::Faithful>(double y);
template double Cbrt<Rounding::Correct>(double y);
}  // namespace method_3²ᴄZ5¹

namespace method_5²Z4¹FMA {

// No overflow or underflow occurs in intermediate computations for
// y ∈ [y₁, y₂].  Steps 2 and 4 respectively involve the 6th powers of q and x,
// approximately the 6th power of ∛y.  The 7th power of ∛y does not overflow for
// y ∈ [y₁, y₂].  Here the range [y₁, y₂] is large enough that a single
// rescaling suffices.
constexpr double y₁ = 0x1p-438;
constexpr double σ₁ = 0x1p-292;
constexpr double σ₁⁻³ = 1 / (σ₁ * σ₁ * σ₁);
constexpr double y₂ = 0x1p438;
constexpr double σ₂ = 0x1p292;
constexpr double σ₂⁻³ = 1 / (σ₂ * σ₂ * σ₂);
static_assert(σ₁⁻³ * y₁ == y₂);
static_assert(σ₂⁻³ * y₂ == y₁);
static_assert(σ₁⁻³ * std::numeric_limits<double>::denorm_min() > y₁);
static_assert(σ₂⁻³ * std::numeric_limits<double>::max() < y₂);

template<Rounding rounding>
double Cbrt(double y) {
  OSACA_FUNCTION_BEGIN(y);
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(masks::sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(masks::sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);

  OSACA_IF (y != y) {
    // The usual logic will produce a qNaN when given a NaN, but will not
    // preserve the payload and will signal overflows (q will be a nonsensical
    // large value, and q³ will overflow).  Further, the rescaling comparisons
    // will signal the invalid operation exception for quiet NaNs (although that
    // would be easy to work around using the unordered compare intrinsics).
    OSACA_RETURN( y + y);
  }

  OSACA_IF(abs_y < y₁) {
    OSACA_IF(abs_y == 0) {
      OSACA_RETURN( y);
    }
    OSACA_RETURN( Cbrt<rounding>(y * σ₁⁻³) * σ₁);
  } OSACA_ELSE_IF(abs_y > y₂) {
    OSACA_IF (abs_y == std::numeric_limits<double>::infinity()) {
      OSACA_RETURN(y);
    }
    OSACA_RETURN( Cbrt<rounding>(y * σ₂⁻³) * σ₂);
  }

  // Step 1.  The constant Γᴋ minimizes the error of q, as in [KB01], rather
  // than that of ξ.  This does not matter all that much here.
  constexpr double Γᴋ = 0x0.19D9'06CB'2868'81F4'88FD'38DF'E7F6'98DD'Bp0;
  constexpr std::uint64_t C = 0x2A9F'7625'3119'D328;
  // The fixed-point number C = ⟦(2×1023 − Γᴋ) / 3⟧ is not representable as a
  // double; frac C + 1 = ⟦2 − Γᴋ / 3⟧ is, and has the same last place.
  constexpr double frac_C_plus_1 = (C & 0x000F'FFFF'FFFF'FFFF) * 0x1p-52 + 1;
  // By sheer luck it happens that ⟦2 − Γᴋ / 3⟧ = ⟦2 − ⟦⟦Γᴋ⟧ / 3⟧⟧.
  static_assert(frac_C_plus_1 == 2 - Γᴋ / 3);
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = _mm_cvtsd_f64(_mm_castsi128_pd(_mm_cvtsi64_si128(Q)));

  // Step 2, the generalized Lagny quadratic irrational method of order 5, with
  // the evaluation strategy described in appendix D.
  double const y² = y * y;
  double const q² = q * q;
  double const q³ = q² * q;
  double const d =
      q / FusedMultiplySubtract(0x1.4A7E9CB8A3491p2 * q, q²,
                                /*-*/ 0x1.08654A2D4F6DBp-1 * abs_y);
  double const ξ = FusedMultiplyAdd(
      d, Sqrt(FusedMultiplySubtract(
                  q³, FusedNegatedMultiplyAdd(q², q, 118.0 / 5 * abs_y),
                  /*-*/ y²)),
      /*+*/ FusedMultiplySubtract(q², q, /*-*/ abs_y) *
            0x1.4A7E9CB8A3491p0 * d);

  // Step 3.
  double const x = _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_sd(ξ), masks::round_toward_zero_26_bits));

  // Step 4, the Lagny–Schröder rational method of order 4.
  double const x² = x * x;  // Exact.
  double const x³ = x² * x;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator = x_sign_y * FusedMultiplySubtract(x², x, abs_y);
  double const denominator =
      FusedMultiplyAdd(x³, FusedMultiplyAdd(10 * x, x², 16 * abs_y), y²);
  double const Δ₁ = FusedMultiplyAdd(6 * x, x², 3 * abs_y);
  double const Δ₂ = numerator / denominator;
  double const r₀ = FusedNegatedMultiplyAdd(Δ₁, Δ₂, x_sign_y);
  double const r₁ = FusedNegatedMultiplyAdd(Δ₁, Δ₂, x_sign_y - r₀);

  double const r̃ = r₀ + 2 * r₁;
  OSACA_IF(rounding == Rounding::Correct &&
      CorrectionPossiblyNeeded(r₀, r₁, r̃, /*τ=*/0x1.E45E16EF5480Fp-76)) {
    OSACA_RETURN(_mm_cvtsd_f64(_mm_or_pd(
        _mm_set_sd(CorrectLastBit(abs_y, std::abs(r₀), std::abs(r̃))), sign)));
  }
  OSACA_RETURN(r₀);
}
template double Cbrt<Rounding::Faithful>(double y);
template double Cbrt<Rounding::Correct>(double y);
}  // namespace method_5²Z4¹FMA

double Cbrt(double const y) {
  return UseHardwareFMA ? method_5²Z4¹FMA::Cbrt<Rounding::Correct>(y)
                        : method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(y);
}

bool CbrtOneBit(double const y, double const a, double const b) {
  double const b² = b * b;
  double const b³ = b² * b;
  DoublePrecision<double> const a² = TwoProduct(a, a);
  auto const& [a²₀, a²₁] = a²;
  DoublePrecision<double> const a³₀ = TwoProduct(a²₀, a);
  DoublePrecision<double> const minus_a³₁ = TwoProduct(a²₁, -a);
  auto const& [a³₀₀, a³₀₁] = a³₀;
  // In cbrt.pdf, where we are specifically considering the computation of the
  // 54th bit, ρ is referred to as ρ₅₃, and ρ_next as ρ₅₄ˌ₁.
  // ρ = y - a³ = y - a³₀ - a³₁ = y - a³₀₀ - a³₀₁ - a³₁;
  double const ρ₀ = y - a³₀₀;  // Exact.
  // ρ = ρ₀ - a³₀₁ - a³₁;
  std::array<double, 4> const ρ = PriestNievergeltNormalize(
      NievergeltQuadruplyCompensatedStep(TwoDifference(ρ₀, a³₀₁), minus_a³₁));
  DCHECK_EQ(ρ[3], 0);
  std::array<double, 3> ρ_next{ρ[0], ρ[1], ρ[2]};
  double const a²₀b = a²₀ * b;
  double const a²₁b = a²₁ * b;
  double const ab² = a * b²;
  for (double rhs : {2 * a²₀b, a²₀b, 2 * a²₁b, a²₁b,  // 3 a²b
                     2 * ab², ab²,                    // 3 ab²
                     b³}) {
    auto const ρ = PriestNievergeltNormalize(NievergeltQuadruplyCompensatedStep(
        TwoSum(ρ_next[0], ρ_next[1]), TwoDifference(ρ_next[2], rhs)));
    DCHECK_EQ(ρ[3], 0);
    ρ_next = {ρ[0], ρ[1], ρ[2]};
  }
  bool const ρ_next_positive =
      ρ_next[0] > 0 || (ρ_next[0] == 0 && ρ_next[1] > 0) ||
      (ρ_next[0] == 0 && ρ_next[1] == 0 && ρ_next[2] >= 0);
  return ρ_next_positive;
}

}  // namespace internal
}  // namespace _cbrt
}  // namespace numerics
}  // namespace principia
