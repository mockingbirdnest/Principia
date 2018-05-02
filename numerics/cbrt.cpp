
#include "numerics/cbrt.hpp"

#include <pmmintrin.h>

#include <cstdint>
#include <limits>

namespace principia {
namespace numerics {

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
