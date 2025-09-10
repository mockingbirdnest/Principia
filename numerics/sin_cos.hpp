#pragma once

#include <functional>

#include "base/macros.hpp"  // 🧙 For FORCE_INLINE.
#include "numerics/fma.hpp"

namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

using namespace principia::numerics::_fma;

using SlowPathCallback = std::function<void(double θ)>;

template<FMAPresence fma_presence>
double __cdecl Sin(double x);
extern template double __cdecl Sin<FMAPresence::Absent>(double x);
extern template double __cdecl Sin<FMAPresence::Present>(double x);

template<FMAPresence fma_presence>
double __cdecl Cos(double x);
extern template double __cdecl Cos<FMAPresence::Absent>(double x);
extern template double __cdecl Cos<FMAPresence::Present>(double x);

template<typename T>
struct SC {
  T sin;
  T cos;
};

template<FMAPresence fma_presence>
SC<double> __cdecl SinCos(double x);
extern template SC<double> __cdecl SinCos<FMAPresence::Absent>(double x);
extern template SC<double> __cdecl SinCos<FMAPresence::Present>(double x);

void SetSlowPathsCallbacks(SlowPathCallback sin_cb,
                           SlowPathCallback cos_cb);

}  // namespace internal

using internal::Cos;
using internal::SetSlowPathsCallbacks;
using internal::Sin;
using internal::SinCos;

}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
