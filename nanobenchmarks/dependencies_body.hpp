#pragma once

#include "nanobenchmarks/dependencies.hpp"

#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_MSVC.
#include "quantities/si.hpp"

#if PRINCIPIA_COMPILER_MSVC
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

namespace principia {
namespace nanobenchmarks {
namespace _dependencies {
namespace internal {

using namespace principia::quantities::_si;

inline double Dependencies<double, double>::ProduceArgument(double const x) {
  return x;
}

inline double Dependencies<double, double>::ConsumeValue(double const value) {
  return value;
}

template<typename Frame>
static Instant Dependencies<Displacement<Frame>, Instant>::ProduceArgument(
    double const x) {
  static Instant const t0;
  return t0 + x * Second;
}

template<typename Frame>
static double Dependencies<Displacement<Frame>, Instant>::ConsumeValue(
    Displacement<Frame> const& value) {
  auto const& coordinates = value.coordinates();
  return _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_pd(coordinates.x / Metre, coordinates.y / Metre),
                 _mm_set_sd(coordinates.z / Metre)));
}

}  // namespace internal

using internal::Dependencies;

}  // namespace _dependencies
}  // namespace nanobenchmarks
}  // namespace principia
