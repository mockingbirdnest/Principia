#pragma once

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
  // The alignment must be sufficiently large to avoid sanitizer complaints.
  return reinterpret_cast<T>(0xBADC0FFEEBADDEC0);
#endif
}

}  // namespace internal
}  // namespace _make_not_null
}  // namespace testing_utilities
}  // namespace principia
