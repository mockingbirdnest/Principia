#pragma once

#include "boost/multiprecision/cpp_bin_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"

namespace principia {
namespace functions {
namespace _multiprecision {
namespace internal {

using namespace boost::multiprecision;

cpp_bin_float_50 Sin(cpp_rational const& α);

}  // namespace internal

using internal::Sin;

}  // namespace _multiprecision
}  // namespace functions
}  // namespace principia
