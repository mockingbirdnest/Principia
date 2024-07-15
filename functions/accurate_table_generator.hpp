#pragma once

#include <functional>
#include <vector>

#include "absl/status/statusor.h"
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

template<typename ArgValue, int degree>
using AccuratePolynomial =
    PolynomialInMonomialBasis<ArgValue, ArgValue, degree>;

template<std::int64_t zeroes>
cpp_rational GalExhaustiveSearch(std::vector<AccurateFunction> const& functions,
                                 cpp_rational const& starting_argument);

template<std::int64_t zeroes>
std::vector<cpp_rational> GalExhaustiveMultisearch(
    std::vector<AccurateFunction> const& functions,
    std::vector<cpp_rational> const& starting_arguments);

// Searches in an interval of radius |T / N| centered on |starting_argument|.
// The |polynomials| must be the degree-2 Taylor approximations of the
// |functions| and the |rests| must be upper bounds on the rest of the Taylor
// series. The argument and function values must be within [1/2, 1[.
template<std::int64_t zeroes>
absl::StatusOr<cpp_rational> StehléZimmermannSimultaneousSearch(
    std::array<AccurateFunction, 2> const& functions,
    std::array<AccuratePolynomial<cpp_rational, 2>, 2> const& polynomials,
    std::array<AccurateFunction, 2> const& rests,
    cpp_rational const& starting_argument,
    std::int64_t N,
    std::int64_t T);

// Performs a search around |starting_argument| to find a solution,
// automatically adjusting the interval over which the search happens.  The
// argument and function values must be nonzero.
template<std::int64_t zeroes>
absl::StatusOr<cpp_rational> StehléZimmermannSimultaneousFullSearch(
    std::array<AccurateFunction, 2> const& functions,
    std::array<AccuratePolynomial<cpp_rational, 2>, 2> const& polynomials,
    std::array<AccurateFunction, 2> const& rests,
    cpp_rational const& starting_argument);

template<std::int64_t zeroes>
std::vector<absl::StatusOr<cpp_rational>>
StehléZimmermannSimultaneousMultisearch(
    std::array<AccurateFunction, 2> const& functions,
    std::vector<std::array<AccuratePolynomial<cpp_rational, 2>, 2>> const&
        polynomials,
    std::vector<std::array<AccurateFunction, 2>> const& rests,
    std::vector<cpp_rational> const& starting_arguments);

}  // namespace internal

using internal::AccurateFunction;
using internal::AccuratePolynomial;
using internal::GalExhaustiveMultisearch;
using internal::GalExhaustiveSearch;
using internal::StehléZimmermannSimultaneousFullSearch;
using internal::StehléZimmermannSimultaneousMultisearch;
using internal::StehléZimmermannSimultaneousSearch;

}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia

#include "functions/accurate_table_generator_body.hpp"
