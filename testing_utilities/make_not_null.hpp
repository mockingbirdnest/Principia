#pragma once

#include "base/not_null.hpp"

namespace principia {

using base::not_null;

namespace testing_utilities {

// Returns a bogus not-null value of the right type.  Dereferencing it is UB.
template<typename T>
not_null<T> make_not_null();

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/make_not_null_body.hpp"
