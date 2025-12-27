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

//template<typename Frame>
//Displacement<Frame> Dependencies<Displacement<Frame>, Instant>::Run(
//    Instant const argument) {
//  double const x = (argument - t0_) / Second;
//  return Displacement<Frame>({x * Metre, x * Metre, x * Metre});
//}
#pragma comment( \
    linker,      \
    "/alternatename:?Run@?$Dependencies@V?$Multivector@V?$Quantity@U?$Dimensions@$00$0A@$0A@$0A@$0A@$0A@$0A@$0A@@internal@_dimensions@quantities@principia@@@internal@_quantities@quantities@principia@@U?$Frame@W4Frame_TestTag@serialization@principia@@$0A@$00$00@2_frame@geometry@5@$00@internal@_grassmann@geometry@principia@@V?$Point@V?$Quantity@U?$Dimensions@$0A@$0A@$00$0A@$0A@$0A@$0A@$0A@@internal@_dimensions@quantities@principia@@@internal@_quantities@quantities@principia@@@2_point@45@@internal@_dependencies@nanobenchmarks@principia@@SQ?AV?$Multivector@V?$Quantity@U?$Dimensions@$00$0A@$0A@$0A@$0A@$0A@$0A@$0A@@internal@_dimensions@quantities@principia@@@internal@_quantities@quantities@principia@@U?$Frame@W4Frame_TestTag@serialization@principia@@$0A@$00$00@2_frame@geometry@5@$00@2_grassmann@geometry@5@V?$Point@V?$Quantity@U?$Dimensions@$0A@$0A@$00$0A@$0A@$0A@$0A@$0A@@internal@_dimensions@quantities@principia@@@internal@_quantities@quantities@principia@@@2_point@85@@Z=foofoofoo")

template<typename Frame>
double Dependencies<Displacement<Frame>, Instant>::ConsumeValue(
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
