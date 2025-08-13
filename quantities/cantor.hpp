#pragma once

#include <concepts>

#include "base/traits.hpp"
#include "boost/multiprecision/number.hpp"

namespace principia {
namespace quantities {
namespace _cantor {
namespace internal {

using namespace boost::multiprecision;
using namespace principia::base::_traits;

// The `cpp_` concepts should be used sparingly, and only in places where the
// Boost multiprecision API differs from our or from the standard C++ API.  At
// any rate, the multiprecision traits should never be used directly.

template<typename T>
concept cpp_number = (is_number<T>::value || is_number_expression<T>::value) &&
                     number_category<T>::value != number_kind_unknown;

template<typename T>
concept cpp_bin_float =
    cpp_number<T> && number_category<T>::value == number_kind_floating_point;

template<typename T>
concept cpp_int =
    cpp_number<T> && number_category<T>::value == number_kind_integer;

template<typename T>
concept cpp_rational =
    cpp_number<T> && number_category<T>::value == number_kind_rational;


template<typename T>
concept discrete = std::integral<T> || cpp_int<T>;

template<typename T>
concept countable = discrete<T> || cpp_rational<T>;

template<typename T>
concept continuum = std::floating_point<T> || cpp_bin_float<T>;

}  // namespace internal

using internal::continuum;
using internal::countable;
using internal::cpp_bin_float;
using internal::cpp_int;
using internal::cpp_number;
using internal::cpp_rational;
using internal::discrete;

}  // namespace _cantor
}  // namespace quantities
}  // namespace principia
