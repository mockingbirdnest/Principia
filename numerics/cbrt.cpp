
#include "numerics/cbrt.hpp"

#include <pmmintrin.h>

#include <array>
#include <cstdint>
#include <limits>
#include <utility>

#include "numerics/double_precision.hpp"

namespace principia {
namespace numerics {

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
// Nievergelt’s procedure, but with our indexing it is is k + 1 at the return
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

bool CorrectionPossiblyNeeded(double const r₀, double const r₁, double const τ) {
  double const r̃ = r₀ + 2 * r₁;
  // TODO(egg): Do we need the right-hand side of the conjunction here?
  return std::abs(0.5 * (r̃ - r₀) - r₁) <= τ * r₀ && r̃ != r₀;
}

double CorrectLastBit(double const y, double const r₀, double const r₁, double const τ) {
  double const r̃ = r₀ + 2 * r₁;
  if (std::abs(0.5 * (r̃ - r₀) - r₁) > τ * r₀ || r̃ == r₀) {
    return r₀;
  }
  // TODO(egg): Handle negative y.
  CHECK_GT(y, 0);
  double const a = std::min(r₀, r̃);
  double const b = 0.5 * (std::max(r₀, r̃) - a);
  double const b² = b * b;
  double const b³ = b² * b;
  DoublePrecision<double> const a² = TwoProduct(a, a);
  auto const& [a²₀, a²₁] = a²;
  DoublePrecision<double> const a³₀ = TwoProduct(a²₀, a);
  DoublePrecision<double> minus_a³₁ = TwoProduct(a²₁, -a);
  auto const& [a³₀₀, a³₀₁] = a³₀;
  // ρ₅₃ = y - a³ = y - a³₀ - a³₁ = y - a³₀₀ - a³₀₁ - a³₁;
  double const ρ₀ = y - a³₀₀;  // Exact.
  // ρ₅₃ = ρ₀ - a³₀₁ - a³₁;
  std::array<double, 4> const ρ₅₃ = PriestNievergeltNormalize(
      NievergeltQuadruplyCompensatedStep(TwoDifference(ρ₀, a³₀₁), minus_a³₁));
  CHECK_EQ(ρ₅₃[3], 0);
  std::array<double, 3> ρ₅₄{ρ₅₃[0], ρ₅₃[1], ρ₅₃[2]};
  for (double rhs : {2 * a²₀ * b, a²₀ * b, 2 * a²₁ * b, a²₁ * b,  // 3 a²b
                     2 * a * b², a * b²,                          // 3 ab²
                     b³}) {
    auto const ρ = PriestNievergeltNormalize(NievergeltQuadruplyCompensatedStep(
        TwoSum(ρ₅₄[0], ρ₅₄[1]), TwoDifference(ρ₅₄[2], rhs)));
    CHECK_EQ(ρ[3], 0);
    ρ₅₄ = {ρ[0], ρ[1], ρ[2]};
  }
  bool const ρ₅₄_positive = ρ₅₄[0] > 0 || (ρ₅₄[0] == 0 && ρ₅₄[1] > 0) ||
                            (ρ₅₄[0] == 0 && ρ₅₄[1] == 0 && ρ₅₄[2] >= 0);
  return ρ₅₄_positive ? std::max(r₀, r̃) : a;
}

constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sign_exponent_and_sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// No overflow or underflow occurs in intermediate computations for
// y ∈ [y₁, y₂].
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double y₁ = 0x1p-225;
constexpr double σ₁ = 0x1p-154;
constexpr double σ₁⁻³ = 1 / (σ₁ * σ₁ * σ₁);
constexpr double y₂ = 0x1p237;
constexpr double σ₂ = 0x1p154;
constexpr double σ₂⁻³ = 1 / (σ₂ * σ₂ * σ₂);
static_assert(σ₁⁻³ * y₁ == y₂, "Incorrect σ₁");
static_assert(σ₂⁻³ * y₂ == y₁, "Incorrect σ₂");
double Cbrt(double const y) {
  __m128d const y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, y_0);
  __m128d const abs_y_0 = _mm_andnot_pd(sign_bit, y_0);
  double const abs_y = _mm_cvtsd_f64(abs_y_0);
  if (y != y) {
    // The usual logic will produce a qNaN when given a NaN, but will not
    // preserve the payload and will signal overflows (q will be a nonsensical
    // large value, and q³ will overflow).  Further, the rescaling comparisons
    // will signal the invalid operation exception for quiet NaNs (although that
    // would be easy to work around using the unordered compare intrinsics).
    return y + y;
  }
  // TODO(egg): we take the absolute value two or three times when going through
  // the rescaling paths; consider having a cbrt_positive function, or a
  // cbrt_positive_unscaled function and four rescaling paths.
  if (abs_y < y₁) {
    if (abs_y == 0) {
      return y;
    }
    return Cbrt(y * σ₁⁻³) * σ₁;
  } else if (abs_y > y₂) {
    if (abs_y == std::numeric_limits<double>::infinity()) {
      return y;
    }
    return Cbrt(y * σ₂⁻³) * σ₂;
  }
  // Approximate ∛y with an error below 3,2 %.  The value of C is chosen to
  // minimize the maximal error of ξ as an approximation of ∛y, ignoring
  // rounding.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(abs_y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = _mm_cvtsd_f64(_mm_castsi128_pd(_mm_cvtsi64_si128(Q)));
  double const q³ = q * q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ = q - (q³ - abs_y) * q / (2 * q³ + abs_y);
  double const x = _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_sd(ξ), sign_exponent_and_sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  return x_sign_y - numerator / denominator;
}

}  // namespace numerics
}  // namespace principia
