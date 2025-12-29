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
Instant Dependencies<Displacement<Frame>, Instant>::ProduceArgument(
    double const x) {
  return t0_ + x * Second;
}

template<typename Frame>
double Dependencies<Displacement<Frame>, Instant>::ConsumeValue(
    Displacement<Frame> const& value) {
  auto const& coordinates = value.coordinates();
  return _mm_cvtsd_f64(
      _mm_and_pd(_mm_set_pd(coordinates.x / Metre, coordinates.y / Metre),
                 _mm_set_sd(coordinates.z / Metre)));
}

template<typename Frame>
Instant Dependencies<RelativeDegreesOfFreedom<Frame>, Instant>::ProduceArgument(
    double const x) {
  return t0_ + x * Second;
}

template<typename Frame>
double Dependencies<RelativeDegreesOfFreedom<Frame>, Instant>::ConsumeValue(
    RelativeDegreesOfFreedom<Frame> const& value) {
  auto const& d = value.displacement().coordinates();
  auto const& v = value.velocity().coordinates();
  return _mm_cvtsd_f64(
      _mm_and_pd(_mm_and_pd(_mm_set_pd(d.x / Metre, d.y / Metre),
                            _mm_set_pd(d.z / Metre, v.x / (Metre / Second))),
                 _mm_set_pd(v.y / (Metre / Second), v.z / (Metre / Second))));
}

}  // namespace internal
}  // namespace _dependencies
}  // namespace nanobenchmarks
}  // namespace principia
