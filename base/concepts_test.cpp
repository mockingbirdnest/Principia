#include "base/concepts.hpp"

#include "boost/multiprecision/cpp_bin_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using namespace boost::multiprecision;
using namespace principia::base::_concepts;

TEST(Concepts, BoostCpp) {
  static_assert(boost_cpp_number<cpp_bin_float_50>);
  static_assert(boost_cpp_number<cpp_int>);
  static_assert(boost_cpp_number<cpp_rational>);
  static_assert(!boost_cpp_number<double>);

  static_assert(boost_cpp_bin_float<cpp_bin_float_50>);
  static_assert(boost_cpp_bin_float<cpp_bin_float_100>);
  static_assert(!boost_cpp_bin_float<cpp_int>);
  static_assert(!boost_cpp_bin_float<cpp_rational>);
  static_assert(!boost_cpp_bin_float<double>);
  static_assert(!boost_cpp_bin_float<int>);

  static_assert(boost_cpp_int<cpp_int>);
  static_assert(!boost_cpp_int<cpp_rational>);

  static_assert(!boost_cpp_rational<cpp_int>);
  static_assert(boost_cpp_rational<cpp_rational>);
}

}  // namespace base
}  // namespace principia
