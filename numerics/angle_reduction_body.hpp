#pragma once

#include "numerics/angle_reduction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "base/macros.hpp"  // 🧙 For PRINCIPIA_USE_SSE3_INTRINSICS.
#include "numerics/elementary_functions.hpp"
#include "numerics/payne_hanek.mathematica.h"
#include "numerics/sin_cos.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_payne_hanek;
using namespace principia::numerics::_sin_cos;
using namespace principia::quantities::_si;

inline constexpr DoublePrecision<double> π_over_2 = []() {
  DoublePrecision<double> result;
  result.value = 0x1.921FB54442D18p0;
  result.error = 0x1.1A62633145C07p-54;
  return result;
}();

template<std::int64_t precision, typename Angle>
void PayneHanekReduction(Angle const& x,
                         DoublePrecision<Angle>& x_reduced,
                         std::int64_t& quadrant) {
  if (x != x) {
    x_reduced = DoublePrecision<Angle>(x);
    quadrant = 0;
    return;
  }

  // This implementation follows [Mul97, section 8.4].  For our purposes it
  // would be a bit more direct to perform a multiplication with 2 / π, but that
  // would require fiddling with the bit indices all over the place.  It's
  // simpler to just follow the book and adjust the result at the end.
  static_assert(std::numeric_limits<double>::radix == 2);
  constexpr std::int64_t p = precision;
  constexpr std::int64_t n = std::numeric_limits<double>::digits;

  int e;
  double X = std::frexp(x / si::Unit<Angle>, &e);
  // Correct the mantissa and exponent to match [Mul97, p. 155].
  e--;
  X = std::scalbn(X, n);
  if (e < -1) {
    DCHECK_LT(Abs(x), 0.5 * si::Unit<Angle>);
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
  // This is not in [Mul97].  The bit numbered 2 - e can only contribute a value
  // sligtly smaller than 8 to `h`.  However, if this bit is 1 and some bits
  // with lower number of the same chunk are 1, then the contribution to `h` can
  // exceed 8.  Therefore, bit 2 - e must be taken into account when computing
  // `h mod 8`.
  std::int64_t const medium_mod_8 = bit_to_index(2 - e);

  // The bits of the chunks must be scaled up by 2^(n + e + 1 + p) to make
  // `Medium(e, p)` an integer.  The computation of `h` then involves a
  // multiplication by 2^(-2n -p).  All together, this results in the scaling
  // below for the chunk.
  double const scale = std::scalbn(1.0, e - n + 1);

  // The most significant chunk has extra bits that we don't need because they
  // are part of `Left(e, p)`.  Drop them: the first bit that is preserved here
  // is the one numbered n - e + 1.
  double const medium_first_mod = std::scalbn(1.0, n - e + 2);

  DoublePrecision<double> h;
  for (std::int64_t i = medium_last; i >= medium_first; --i) {
    double chunk = PayneHanekChunks[i];
    if (i == medium_first) {
      chunk = std::remainder(chunk, medium_first_mod);
    }
    double const schunk = scale * chunk;
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
      if (i < medium_mod_8) {
        Xl_schunk = std::remainder(Xl_schunk, 8.0);
      } else {
        DCHECK_LE(Xl_schunk, 8.0);
      }
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
  x_reduced = (h - round_h) * (π_over_2 * si::Unit<Angle>);
}

template<>
inline Angle ReduceAngle<-π / 2, π / 2>(Angle const& θ) {
  double θ_reduced;
  std::int64_t quadrant;
  Reduce(θ / Radian, θ_reduced, quadrant);
  if (quadrant == 1 || quadrant == 3) {
    if (θ_reduced < 0.0) {
      θ_reduced += π / 2;
    } else {
      θ_reduced -= π / 2;
    }
  }
  return θ_reduced * Radian;
}

template<>
inline Angle ReduceAngle<-π, π>(Angle const& θ) {
  double θ_reduced;
  std::int64_t quadrant;
  Reduce(θ / Radian, θ_reduced, quadrant);
  if (quadrant == 1) {
    θ_reduced += π / 2;
  } else if (quadrant == 3) {
    θ_reduced -= π / 2;
  } else if (quadrant == 2) {
    if (θ_reduced < 0.0) {
      θ_reduced += π;
    } else {
      θ_reduced -= π;
    }
  }
  return θ_reduced * Radian;
}

template<>
inline Angle ReduceAngle<0.0, 2 * π>(Angle const& θ) {
  double θ_reduced;
  std::int64_t quadrant;
  Reduce(θ / Radian, θ_reduced, quadrant);
  if (quadrant == 1) {
    θ_reduced += π / 2;
  } else if (quadrant == 3) {
    θ_reduced += 3 * π / 2;
  } else if (quadrant == 2) {
    θ_reduced += π;
  } else {
    // `quadrant == 0`.
    if (θ_reduced < 0.0) {
      θ_reduced += 2 * π;
    }
  }
  return θ_reduced * Radian;
}

// θ = fractional_part + integer_part * π where fractional_part is in
// [-π/2, π/2].
template<>
inline void ReduceAngle<-π / 2, π / 2>(Angle const& θ,
                                       Angle& fractional_part,
                                       double& integer_part) {
  fractional_part = ReduceAngle<-π / 2, π / 2>(θ);
  integer_part = std::round((θ - fractional_part) / (π * Radian));
}

// θ = fractional_part + integer_part * 2 * π where fractional_part is in
// [-π, π].
template<>
inline void ReduceAngle<-π, π>(Angle const& θ,
                               Angle& fractional_part,
                               double& integer_part) {
  fractional_part = ReduceAngle<-π, π>(θ);
  integer_part = std::round((θ - fractional_part) / (2 * π * Radian));
}

// θ = fractional_part + integer_part * 2 * π where fractional_part is in
// [0, 2 * π].
template<>
inline void ReduceAngle<0.0, 2 * π>(Angle const& θ,
                                    Angle& fractional_part,
                                    double& integer_part) {
  fractional_part = ReduceAngle<0.0, 2 * π>(θ);
  integer_part = std::round((θ - fractional_part) / (2 * π * Radian));
}

}  // namespace internal
}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia
