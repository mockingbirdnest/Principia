#pragma once

#include "numerics/angle_reduction.hpp"

#include <cmath>
#include <limits>

#include "base/macros.hpp"  // 🧙 For PRINCIPIA_USE_SSE3_INTRINSICS.
#include "numerics/double_precision.hpp"
#include "numerics/payne_hanek.mathematica.h"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::numerics::_double_precision;
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
  // This implementation follows [Mul97, section 8.4].
  static_assert(std::numeric_limits<double>::radix == 2);
  constexpr std::int64_t p = precision;
  constexpr std::int64_t n = std::numeric_limits<double>::digits;

  int e;
  double X = std::frexp(x / Radian, &e);
  // Correct the mantissa and exponent to match [Mul97, p. 155].
  e--;
  X = std::scalbn(X, n);
  CHECK_GE(e, -1);

  // Split `X` so that the multiplications below are exact.
  double Xh;
  double Xl;
  VeltkampSplitting<PayneHanekBitsPerChunk>(X, Xh, Xl);

  // Converts the bit numbering of [Mul97, p. 155] into our index in
  // `PayneHanekChunks`.
  auto bit_to_index = [](std::int64_t const b) {
    return -b / PayneHanekBitsPerChunk;
  };

  // `Medium(e, p)` is made of the chunks with indices [medium_first,
  // medium_last].
  std::int64_t const medium_first = bit_to_index(n - e + 1);
  std::int64_t const medium_last = bit_to_index(-n - e - 1 - p);

  DoublePrecision<double> h;
  double h_quadrant = 0.0;
  for (std::int64_t i = medium_last; i >= medium_first; --i) {
    double const chunk = std::scalbn(PayneHanekChunks[i], n + e + 1 + p);
    double const Xl_chunk = Xl * chunk;
    double const Xh_chunk = Xh * chunk;
    h_quadrant += std::remainder(Xl_chunk, 4.0) + std::remainder(Xh_chunk, 4.0);
    h += Xl_chunk;
    h += Xh_chunk;
  }

  h -= std::trunc(h.value);
  h -= std::trunc(h.value);
  quadrant = static_cast<std::int64_t>(h_quadrant) % 4;
  x_reduced = Scale(Radian, h);
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
