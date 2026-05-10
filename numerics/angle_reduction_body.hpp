#pragma once

#include "numerics/angle_reduction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "base/macros.hpp"  // 🧙 For PRINCIPIA_USE_SSE3_INTRINSICS.
#include "numerics/double_precision.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/payne_hanek.mathematica.h"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_payne_hanek;
using namespace principia::quantities::_si;

// A somewhat arbitrary value above which we fail argument reduction.
// TODO(phl): We *really* need a proper angle reduction.
double const reduction_threshold =
    static_cast<double>(std::numeric_limits<std::int64_t>::max() / 2);

template<>
inline constexpr Angle one_π<Angle> = π * Radian;

template<>
inline constexpr DoublePrecision<Angle> one_π<DoublePrecision<Angle>> = []() {
  DoublePrecision<Angle> result;
  result.value = 0x1.921FB54442D18p1 * Radian;
  result.error = 0x1.1A62633145C07p-53 * Radian;
  return result;
}();

template<>
inline constexpr Angle two_π<Angle> = 2 * π * Radian;

template<>
inline constexpr DoublePrecision<Angle> two_π<DoublePrecision<Angle>> = []() {
  DoublePrecision<Angle> result;
  result.value = 0x1.921FB54442D18p2 * Radian;
  result.error = 0x1.1A62633145C07p-52 * Radian;
  return result;
}();

inline constexpr DoublePrecision<Angle> π_over_2 = []() {
  DoublePrecision<Angle> result;
  result.value = 0x1.921FB54442D18p0 * Radian;
  result.error = 0x1.1A62633145C07p-54 * Radian;
  return result;
}();

template<typename T>
std::int64_t StaticCastToInt64(T const& t);

template<>
inline std::int64_t StaticCastToInt64<double>(double const& t) {
  return static_cast<std::int64_t>(t);
}

template<>
inline std::int64_t StaticCastToInt64<DoublePrecision<double>>(
    DoublePrecision<double> const& t) {
  return static_cast<std::int64_t>(t.value + t.error);
}

template<typename Angle,
         DoubleWrapper fractional_part_lower_bound,
         DoubleWrapper fractional_part_upper_bound>
class AngleReduction;

// TODO(phl): This is extremely imprecise near large multiples of π.  Use a
// better algorithm (Payne-Hanek?).
template<>
class AngleReduction<Angle, -π / 2, π / 2> {
 public:
  // Argument reduction: angle = fractional_part + integer_part * π where
  // fractional_part is in [-π/2, π/2].
  static bool Reduce(Angle const& θ,
                     Angle& fractional_part,
                     std::int64_t& integer_part) {
    if (θ > reduction_threshold * one_π<Angle> ||
        θ < -reduction_threshold * one_π<Angle>) {
      return false;
    }
    double const θ_in_half_cycles = θ / (π * Radian);
    double reduced_in_half_cycles;
#if PRINCIPIA_USE_SSE3_INTRINSICS
    auto const& x = θ_in_half_cycles;
    __m128d const x_128d = _mm_set_sd(x);
    integer_part = _mm_cvtsd_si64(x_128d);
    reduced_in_half_cycles = _mm_cvtsd_f64(
        _mm_sub_sd(x_128d, _mm_cvtsi64_sd(__m128d{}, integer_part)));
#else
    integer_part = static_cast<std::int64_t>(std::nearbyint(θ_in_half_cycles));
    reduced_in_half_cycles = θ_in_half_cycles - integer_part;
#endif
    fractional_part = reduced_in_half_cycles * π * Radian;
    return true;
  }
};

template<typename Angle>
class AngleReduction<Angle, -π, π> {
 public:
  static bool Reduce(Angle const& θ,
                     Angle& fractional_part,
                     std::int64_t& integer_part) {
    if (!AngleReduction<Angle, 0.0, 2 * π>::Reduce(
            θ, fractional_part, integer_part)) {
      return false;
    }
    if (fractional_part > one_π<Angle>) {
      fractional_part -= two_π<Angle>;
      ++integer_part;
    }
    return true;
  }
};

template<typename Angle>
class AngleReduction<Angle, 0.0, 2 * π> {
 public:
  static bool Reduce(Angle const& θ,
                     Angle& fractional_part,
                     std::int64_t& integer_part) {
    if (!AngleReduction<Angle, -2 * π, 2 * π>::Reduce(
            θ, fractional_part, integer_part)) {
      return false;
    }
    if (fractional_part < Angle{}) {
      fractional_part += two_π<Angle>;
      --integer_part;
    }
    return true;
  }
};

template<typename Angle>
class AngleReduction<Angle, -2 * π, 2 * π> {
 public:
  static bool Reduce(Angle const& θ,
                     Angle& fractional_part,
                     std::int64_t& integer_part) {
    if (θ > reduction_threshold * one_π<Angle> ||
        θ < -reduction_threshold * one_π<Angle>) {
      return false;
    }
    // This has the same semantics as fmod.
    auto const θ_over_2π = θ / two_π<Angle>;
    integer_part = StaticCastToInt64(θ_over_2π);
    fractional_part =
        θ - two_π<Angle> * static_cast<decltype(θ_over_2π)>(integer_part);
    return true;
  }
};

template<std::int64_t precision>
void PayneHanek(Angle const& x,
                DoublePrecision<Angle>& x_reduced,
                std::int64_t& quadrant) {
  // This implementation follows [Mul97, section 8.4].  For our purposes it
  // would be a bit more direct to perform a multiplication with 2 / π, but that
  // would require fiddling with the bit indices all over the place.  It's
  // simpler to just follow the book and adjust the result at the end.
  static_assert(std::numeric_limits<double>::radix == 2);
  constexpr std::int64_t p = precision;
  constexpr std::int64_t n = std::numeric_limits<double>::digits;

  int e;
  double X = std::frexp(x / Radian, &e);
  // Correct the mantissa and exponent to match [Mul97, p. 155].
  e--;
  X = std::scalbn(X, n);
  if (e < -1) {
    DCHECK_LT(Abs(x), 0.5 * Radian);
    x_reduced = DoublePrecision<Angle>(x);
    quadrant = 0;
    return;
  }

  // Split `X` so that the multiplications below are exact.
  double Xh;
  double Xl;
  VeltkampSplitting<PayneHanekBitsPerChunk>(X, Xh, Xl);

  // Converts the bit numbering of [Mul97, p. 155] into our index in
  // `PayneHanekChunks`.
  auto bit_to_index = [](std::int64_t const b) {
    return std::max(INT64_C(0), -b / PayneHanekBitsPerChunk);
  };

  // `Medium(e, p)` is made of the chunks with indices [medium_first,
  // medium_last].
  std::int64_t const medium_first = bit_to_index(n - e + 1);
  std::int64_t const medium_last = bit_to_index(-n - e - 1 - p);
  // This is not in [Mul97].  The bit numbered 2 - e is the smallest one that
  // can contribute to `h mod 8`.
  std::int64_t const medium_mod_8 = bit_to_index(2 - e);

  DoublePrecision<double> h;
  for (std::int64_t i = medium_last; i >= medium_first; --i) {
    double chunk = PayneHanekChunks[i];
    if (i == medium_first) {
      // The most significant chunk has extra bits that we don't need because
      // they are part of `Left(e, p)`.  Drop them: the first bit that is
      // preserved here is the one numbered n - e + 1.
      chunk = std::remainder(chunk, std::scalbn(1.0, n - e + 2));
    }
    // The bits of the chunks must be scaled up by 2^(n + e + 1 + p) to make
    // `Medium(e, p)` an integer.  The computation of `h` then involves a
    // multiplication by 2^(-2n -p).  All together, this results in the scaling
    // below for the chunk.
    double const schunk = std::scalbn(chunk, e - n + 1);
    // The products are exact by construction of the chunks.
    double Xl_schunk = Xl * schunk;
    double Xh_schunk = Xh * schunk;
    // `h` as defined in [Mul97] may be as large as 2^(n + 3).  Since we have
    // at most 2n bits in our `DoublePrecision` result, this would give us at
    // most n - 3 bits after the fractional point, which is not sufficient for
    // all argument reductions.  We reduce the products that can contribute to
    // `h mod 8` before doing the additions.  There are approximately n / 26
    // chunks that can contribute to `h mod 8`, so that gives us an upper bound
    // of roughly 0.308 n for `h`.
    if (i <= medium_mod_8) {
      Xl_schunk = std::remainder(Xl_schunk, 8.0);
      Xh_schunk = std::remainder(Xh_schunk, 8.0);
    } else {
      DCHECK_LE(Xl_schunk, 8.0);
      DCHECK_LE(Xh_schunk, 8.0);
    }
    h += Xl_schunk;
    h += Xh_schunk;
  }

  // This is where we diverge from [Mul97]: we really want x * (2 / π) with
  // rounding, not x * (4 / π) with truncation.
  h = Scale(0.5, h);

  // Because of the bound on `h` above, the fractional point cannot be in the
  // `error`.
  double const round_h = std::round(h.value);
  quadrant = static_cast<std::int64_t>(round_h) & 0b11;
  x_reduced = (h - round_h) * π_over_2;
}

template<DoubleWrapper fractional_part_lower_bound,
         DoubleWrapper fractional_part_upper_bound,
         typename Angle>
bool ReduceAngle(Angle const& θ,
                 Angle& fractional_part,
                 std::int64_t& integer_part) {
  return AngleReduction<Angle,
                        fractional_part_lower_bound,
                        fractional_part_upper_bound>::Reduce(θ,
                                                             fractional_part,
                                                             integer_part);
}

template<DoubleWrapper fractional_part_lower_bound,
         DoubleWrapper fractional_part_upper_bound,
         typename Angle>
Angle ReduceAngle(Angle const& θ) {
  Angle fractional_part;
  std::int64_t integer_part;
  CHECK((AngleReduction<Angle,
                        fractional_part_lower_bound,
                        fractional_part_upper_bound>::Reduce(θ,
                                                             fractional_part,
                                                             integer_part)));
  return fractional_part;
}

}  // namespace internal
}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia
