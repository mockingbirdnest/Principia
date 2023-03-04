#pragma once

#include "base/macros.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace testing_utilities {
namespace _make_not_null {
namespace internal {

template<typename T>
not_null<T> make_not_null() {
#if ARCH_CPU_32_BITS
  return reinterpret_cast<T>(0xDEADBEEF);
#elif ARCH_CPU_64_BITS
  return reinterpret_cast<T>(0xBADC0FFEE0DDF00D);
#endif
}

}  // namespace internal
}  // namespace _make_not_null
}  // namespace testing_utilities
}  // namespace principia
