#pragma once

#include "base/macros.hpp"  // 🧙 For FORCE_INLINE.

namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

double __cdecl Sin(double x);
double __cdecl Cos(double x);

}  // namespace internal

using internal::Cos;
using internal::Sin;

}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
