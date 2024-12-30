#pragma once

#include "base/macros.hpp"

namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

#define PRINCIPIA_INLINE_SIN_COS 1
#define PRINCIPIA_USE_OSACA_SIN 0
#define PRINCIPIA_USE_OSACA_COS 1

#if PRINCIPIA_INLINE_SIN_COS
FORCE_INLINE(inline)
#endif
double __cdecl Sin(double x);
#if PRINCIPIA_INLINE_SIN_COS
FORCE_INLINE(inline)
#endif
double __cdecl Cos(double x);

}  // namespace internal

using internal::Cos;
using internal::Sin;

}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia

#if PRINCIPIA_INLINE_SIN_COS
#include "numerics/sin_cos.cpp"
#endif
