#pragma once

#include <concepts>

#include "base/traits.hpp"
#include "boost/multiprecision/number.hpp"
#include "quantities/m128d.hpp"

namespace principia {
namespace quantities {
namespace _cantor {
namespace internal {

using namespace boost::multiprecision;
using namespace principia::base::_traits;
using namespace principia::quantities::_m128d;

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

template<typename T>
concept discrete = std::integral<T> || boost_cpp_int<T>;

template<typename T>
concept countable = discrete<T> || boost_cpp_rational<T>;

template<typename T>
concept continuum =
    std::floating_point<T> || std::same_as<T, M128D> || boost_cpp_bin_float<T>;

}  // namespace internal

using internal::boost_cpp_bin_float;
using internal::boost_cpp_int;
using internal::boost_cpp_number;
using internal::boost_cpp_rational;
using internal::continuum;
using internal::countable;
using internal::discrete;

}  // namespace _cantor
}  // namespace quantities
}  // namespace principia
