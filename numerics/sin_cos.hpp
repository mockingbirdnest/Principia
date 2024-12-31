#pragma once

#include "base/macros.hpp"  // 🧙 For FORCE_INLINE.

namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

#define PRINCIPIA_INLINE_SIN_COS 0
#define OSACA_ANALYSED_FUNCTION

#if defined(OSACA_ANALYSED_FUNCTION)
#define PRINCIPIA_USE_OSACA !PRINCIPIA_MACRO_IS_EMPTY(OSACA_ANALYSED_FUNCTION)
#endif

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
