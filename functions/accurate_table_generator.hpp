#pragma once

#include <functional>
#include <vector>

#include "boost/multiprecision/cpp_bin_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"

namespace principia {
namespace functions {
namespace _accurate_table_generator {
namespace internal {

using namespace boost::multiprecision;

using AccurateFunction = std::function<cpp_bin_float_50(cpp_rational const&)>;

template<std::int64_t zeroes>
cpp_rational ExhaustiveSearch(std::vector<AccurateFunction> const& functions,
                              cpp_rational const& start);

}  // namespace internal

using internal::ExhaustiveSearch;

}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia

#include "functions/accurate_table_generator_body.hpp"