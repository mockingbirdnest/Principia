#pragma once

namespace principia {
namespace numerics {
namespace _fast_sin_cos_2π {
namespace internal {

// An implementation of sin and cos that favours performance at the expense of
// accuracy.  It has an absolute error of roughly 6e-7 and takes about 25 CPU
// cycles.  The argument must be in the range of the 64-bit integers.
void FastSinCos2π(double cycles, double& sin, double& cos);

}  // namespace internal

using internal::FastSinCos2π;

}  // namespace _fast_sin_cos_2π
}  // namespace numerics
}  // namespace principia

namespace principia::numerics {
using namespace principia::numerics::_fast_sin_cos_2π;
}  // namespace principia::numerics
