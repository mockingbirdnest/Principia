#include "quantities/cantor.hpp"

#include "boost/multiprecision/cpp_bin_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "gtest/gtest.h"

namespace principia {
namespace quantities {

using namespace boost::multiprecision;
using namespace principia::quantities::_cantor;

TEST(Cantor, CppBinFloat) {
  static_assert(boost_cpp_bin_float<cpp_bin_float_50>);
  static_assert(boost_cpp_bin_float<cpp_bin_float_100>);
  static_assert(!boost_cpp_bin_float<cpp_int>);
  static_assert(!boost_cpp_bin_float<cpp_rational>);
  static_assert(!boost_cpp_bin_float<double>);
  static_assert(!boost_cpp_bin_float<int>);
}

TEST(Cantor, Discrete) {
  static_assert(discrete<int>);
  static_assert(discrete<cpp_int>);
  static_assert(!discrete<cpp_rational>);
  static_assert(!discrete<double>);
  static_assert(!discrete<cpp_bin_float_50>);
}

TEST(Cantor, Countable) {
  static_assert(countable<int>);
  static_assert(countable<cpp_int>);
  static_assert(countable<cpp_rational>);
  static_assert(!countable<double>);
  static_assert(!countable<cpp_bin_float_50>);
}

TEST(Cantor, Continuum) {
  static_assert(continuum<double>);
  static_assert(continuum<cpp_bin_float_50>);
  static_assert(!continuum<int>);
  static_assert(!continuum<cpp_int>);
  static_assert(!continuum<cpp_rational>);
}

}  // namespace quantities
}  // namespace principia
