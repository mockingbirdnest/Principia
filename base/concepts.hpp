#pragma once

#include "boost/multiprecision/number.hpp"

namespace principia {
namespace base {
namespace _concepts {
namespace internal {

using namespace boost::multiprecision;

// True if and only if T has a (possibly templated) static member function named
// ReadFromMessage.
template<typename T>
concept serializable = requires {
  &T::ReadFromMessage;
} || requires {
  &T::template ReadFromMessage<>;
};

// The `boost_cpp_` concepts should be used sparingly, and only in places where
// the Boost multiprecision API differs from ours or from the standard C++ API.
// At any rate, the multiprecision traits should never be used directly.

template<typename T>
concept boost_cpp_number =
    (is_number<T>::value || is_number_expression<T>::value) &&
    number_category<T>::value != number_kind_unknown;

template<typename T>
concept boost_cpp_bin_float =
    boost_cpp_number<T> &&
    number_category<T>::value == number_kind_floating_point;

template<typename T>
concept boost_cpp_int =
    boost_cpp_number<T> && number_category<T>::value == number_kind_integer;

template<typename T>
concept boost_cpp_rational =
    boost_cpp_number<T> && number_category<T>::value == number_kind_rational;

}  // namespace internal

using internal::boost_cpp_bin_float;
using internal::boost_cpp_int;
using internal::boost_cpp_number;
using internal::boost_cpp_rational;
using internal::serializable;

}  // namespace _concepts
}  // namespace base
}  // namespace principia
