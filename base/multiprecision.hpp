#pragma once

#include "boost/multiprecision/cpp_bin_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "boost/multiprecision/number.hpp"

namespace principia {
namespace base {
namespace _multiprecision {
namespace internal {

// The `boost_cpp_` concepts should be used sparingly, and only in places where
// the Boost multiprecision API differs from ours or from the standard C++ API.
// At any rate, the multiprecision traits should never be used directly.

template<typename T>
concept boost_cpp_number =
    (boost::multiprecision::is_number<T>::value ||
     boost::multiprecision::is_number_expression<T>::value) &&
    boost::multiprecision::number_category<T>::value !=
        boost::multiprecision::number_kind_unknown;

template<typename T>
concept boost_cpp_bin_float =
    boost_cpp_number<T> &&
    boost::multiprecision::number_category<T>::value ==
        boost::multiprecision::number_kind_floating_point;

template<typename T>
concept boost_cpp_int =
    boost_cpp_number<T> && boost::multiprecision::number_category<T>::value ==
                               boost::multiprecision::number_kind_integer;

template<typename T>
concept boost_cpp_rational =
    boost_cpp_number<T> && boost::multiprecision::number_category<T>::value ==
                               boost::multiprecision::number_kind_rational;

using cpp_rational =
    boost::multiprecision::number<boost::multiprecision::cpp_rational_backend,
                                  boost::multiprecision::et_off>;

using cpp_int =
    boost::multiprecision::number<boost::multiprecision::cpp_int_backend<>,
                                  boost::multiprecision::et_off>;

template<unsigned digits>
using cpp_bin_float = boost::multiprecision::number<
    boost::multiprecision::backends::cpp_bin_float<digits>,
    boost::multiprecision::et_off>;

using cpp_bin_float_50 = cpp_bin_float<50>;

using cpp_bin_float_100 = cpp_bin_float<100>;

}  // namespace internal

using internal::boost_cpp_bin_float;
using internal::boost_cpp_int;
using internal::boost_cpp_number;
using internal::boost_cpp_rational;
using internal::cpp_bin_float;
using internal::cpp_bin_float_100;
using internal::cpp_bin_float_50;
using internal::cpp_int;
using internal::cpp_rational;

}  // namespace _multiprecision
}  // namespace base
}  // namespace principia
