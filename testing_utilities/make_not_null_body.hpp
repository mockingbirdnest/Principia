#pragma once

#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace testing_utilities {

// Returns a bogus not-null value of the right type.  Dereferencing it is UB.
template<typename T>
not_null<T> make_not_null() {
  return reinterpret_cast<T>(0xDEADBEEF);
}

}  // namespace testing_utilities
}  // namespace principia
