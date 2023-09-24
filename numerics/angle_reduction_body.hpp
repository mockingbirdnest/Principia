#pragma once

#include "numerics/angle_reduction.hpp"

namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::quantities::_si;

// TODO(phl): This is extremely imprecise near large multiples of π.  Use a
// better algorithm (Payne-Hanek?).
void Reduce(Angle const& angle,
            Angle& fractional_part,
            std::int64_t& integer_part) {
  double const angle_in_half_cycles = angle / (π * Radian);
  double reduced_in_half_cycles;
#if PRINCIPIA_USE_SSE3_INTRINSICS
  auto const& x = angle_in_half_cycles;
  __m128d const x_128d = _mm_set_sd(x);
  integer_part = _mm_cvtsd_si64(x_128d);
  reduced_in_half_cycles = _mm_cvtsd_f64(
      _mm_sub_sd(x_128d,
                 _mm_cvtsi64_sd(__m128d{}, integer_part)));
#else
  integer_part = std::nearbyint(angle_in_half_cycles);
  reduced_in_half_cycles = angle_in_half_cycles - integer_part;
#endif
  fractional_part = reduced_in_half_cycles * π * Radian;
}

inline DoublePrecision<Angle> Mod2π(DoublePrecision<Angle> const& θ) {
  static DoublePrecision<Angle> const two_π = []() {
    return QuickTwoSum(0x1.921FB54442D18p2 * Radian,
                       0x1.1A62633145C07p-52 * Radian);
  }();
  auto const θ_over_2π = θ / two_π;
  return θ - two_π * DoublePrecision<double>(static_cast<int>(θ_over_2π.value));
}

}  // namespace internal
}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia
