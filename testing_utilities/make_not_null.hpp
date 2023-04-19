#pragma once

#include "base/not_null.hpp"

namespace principia {
namespace testing_utilities {
namespace _make_not_null {
namespace internal {

using namespace principia::base::_not_null;

// Returns a bogus not-null value of the right type.  Dereferencing it is UB.
template<typename T>
not_null<T> make_not_null();

}  // namespace internal

using internal::make_not_null;

}  // namespace _make_not_null
}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/make_not_null_body.hpp"
