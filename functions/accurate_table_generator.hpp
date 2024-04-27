#pragma once

#include <functional>
#include <vector>

#include "boost/multiprecision/cpp_bin_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"

namespace principia {
namespace functions {
namespace _accurate_table_generator {
namespace internal {

using namespace boost::multiprecision;
using namespace principia::numerics::_polynomial_in_monomial_basis;

using AccurateFunction = std::function<cpp_bin_float_50(cpp_rational const&)>;
#if 1
template<int degree>
using AccuratePolynomial =
    PolynomialInMonomialBasis<cpp_rational, cpp_rational, degree>;
#endif

template<std::int64_t zeroes>
cpp_rational ExhaustiveSearch(std::vector<AccurateFunction> const& functions,
                              cpp_rational const& starting_argument);

template<std::int64_t zeroes>
std::vector<cpp_rational> ExhaustiveMultisearch(
    std::vector<AccurateFunction> const& functions,
    std::vector<cpp_rational> const& starting_arguments);

#if 1
template<std::int64_t zeroes>
cpp_rational SimultaneousBadCaseSearch(
  std::array<AccurateFunction, 2> const& functions,
  std::array<AccuratePolynomial<2>, 2> const& polynomials,
  std::int64_t const M,
  std::int64_t const T);
#endif

}  // namespace internal

using internal::AccuratePolynomial;
using internal::ExhaustiveMultisearch;
using internal::ExhaustiveSearch;
#if 1
using internal::SimultaneousBadCaseSearch;
#endif

}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia

#include "functions/accurate_table_generator_body.hpp"
