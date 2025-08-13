#pragma once

#include <concepts>

#include "base/traits.hpp"
#include "boost/multiprecision/detail/number_base.hpp"
#include "boost/multiprecision/fwd.hpp"

namespace principia {
namespace quantities {
namespace _cantor {
namespace internal {

using namespace boost::multiprecision;
using namespace principia::base::_traits;

// This concept should be used sparingly, and only in places where the Boost
// multiprecision API differs from our or from the standard C++ API.
template<typename T>
concept cpp_bin_float = is_number<T>::value &&
                        number_category<T>::value == number_kind_floating_point;

template<typename T>
concept discrete =
    std::integral<T> || std::same_as<T, boost::multiprecision::cpp_int>;

template<typename T>
concept countable =
    discrete<T> || std::same_as<T, boost::multiprecision::cpp_rational>;

template<typename T>
concept continuum = std::floating_point<T> || cpp_bin_float<T>;

}  // namespace internal

using internal::continuum;
using internal::cpp_bin_float;
using internal::countable;
using internal::discrete;

}  // namespace _cantor
}  // namespace quantities
}  // namespace principia
