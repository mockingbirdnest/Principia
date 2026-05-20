#pragma once

#include <functional>

#include "numerics/fma.hpp"

namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

using namespace principia::numerics::_fma;

using SlowPathCallback = std::function<void(double θ)>;

// It would be inconvenient to move angle reduction to `angle_reduction.hpp`
// because of the entanglement with OSACA.  It would also be unpleasant to
// duplicate the code.  Let's just expose it for use by `angle_reduction.hpp`.
// Note that we don't return a `DoublePrecision` as that would cause circular
// dependencies.
void Reduce(double x, double& x_reduced, std::int64_t& quadrant);

template<FMAPresence fma_presence>
double __cdecl Sin(double x);
extern template double __cdecl Sin<FMAPresence::Absent>(double x);
extern template double __cdecl Sin<FMAPresence::Present>(double x);

template<FMAPresence fma_presence>
double __cdecl Cos(double x);
extern template double __cdecl Cos<FMAPresence::Absent>(double x);
extern template double __cdecl Cos<FMAPresence::Present>(double x);

// Keep this type internal and use structured bindings to extract the parts.
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
using internal::Reduce;
using internal::Sin;
using internal::SinCos;

}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
