#pragma once

#include <functional>
#include <vector>

#include "absl/status/statusor.h"
#include "base/thread_pool.hpp"
#include "boost/multiprecision/cpp_bin_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"

namespace principia {
namespace functions {
namespace _accurate_table_generator {
namespace internal {

using namespace boost::multiprecision;
using namespace principia::base::_thread_pool;
using namespace principia::numerics::_polynomial_in_monomial_basis;

using AccurateFunction = std::function<cpp_bin_float_50(cpp_rational const&)>;

// The use of factories below greatly speeds up the search (in one case, from
// multiple hours to 11 s) without affecting correctness.  It seems that Taylor
// polynomials constructed at the `starting_argument` get transformed into
// polynomials with larger and larger coefficients as we scan more and more
// distant slices.  This in turn makes polynomial composition and lattice
// reduction progressively more expensive.  Locally constructed Taylor
// polynomials behave much better (100'000× speed-ups have been observed for
// some slices).

template<typename ArgValue, int degree>
using AccuratePolynomial =
    PolynomialInMonomialBasis<ArgValue, ArgValue, degree>;
template<typename ArgValue, int degree>
using AccuratePolynomialFactory =
    std::function<AccuratePolynomial<ArgValue, degree>(cpp_rational const&)>;

// The remainders don't need to be extremely precise, so for speed
// they are computed using double.
using ApproximateFunction = std::function<double(cpp_rational const&)>;
using ApproximateFunctionFactory =
    std::function<ApproximateFunction(cpp_rational const&)>;


template<std::int64_t zeroes>
cpp_rational GalExhaustiveSearch(std::vector<AccurateFunction> const& functions,
                                 cpp_rational const& starting_argument);

template<std::int64_t zeroes>
std::vector<cpp_rational> GalExhaustiveMultisearch(
    std::vector<AccurateFunction> const& functions,
    std::vector<cpp_rational> const& starting_arguments);

// Searches in an interval of radius `T / N` centered on `starting_argument`.
// The `polynomials` must be the degree-2 Taylor approximations of the
// `functions` and the `remainders` must be upper bounds on the remainder of the
// Taylor series. The argument and function values must be within [1/2, 1[.
template<std::int64_t zeroes>
absl::StatusOr<cpp_rational> StehléZimmermannSimultaneousSearch(
    std::array<AccurateFunction, 2> const& functions,
    std::array<AccuratePolynomial<cpp_rational, 2>, 2> const& polynomials,
    std::array<ApproximateFunction, 2> const& remainders,
    cpp_rational const& starting_argument,
    std::int64_t N,
    std::int64_t T);

// Performs a search around `starting_argument` to find a solution,
// automatically adjusting the interval over which the search happens.  The
// argument and function values must be nonzero.  If `search_pool` is not null,
// the search may use speculative execution.
template<std::int64_t zeroes>
absl::StatusOr<cpp_rational> StehléZimmermannSimultaneousFullSearch(
    std::array<AccurateFunction, 2> const& functions,
    std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2> const&
        polynomials,
    std::array<ApproximateFunctionFactory, 2> const& remainders,
    cpp_rational const& starting_argument,
    ThreadPool<void>* search_pool = nullptr);

// Same as above, but performs searches in parallel using the corresponding
// `polynomials`, `remainders`, and `starting_arguments`.  Returns the results
// in the same order as the parameters.
template<std::int64_t zeroes>
std::vector<absl::StatusOr<cpp_rational>>
StehléZimmermannSimultaneousMultisearch(
    std::array<AccurateFunction, 2> const& functions,
    std::vector<std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2>>
        const& polynomials,
    std::vector<std::array<ApproximateFunctionFactory, 2>> const& remainders,
    std::vector<cpp_rational> const& starting_arguments);

// Same as above, but instead of accumulating all the results and returning them
// in a vector, it runs `callback` each time a computation is complete.  The
// `index` indicates to which parameters the result corresponds.
template<std::int64_t zeroes>
void StehléZimmermannSimultaneousStreamingMultisearch(
    std::array<AccurateFunction, 2> const& functions,
    std::vector<std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2>>
        const& polynomials,
    std::vector<std::array<ApproximateFunctionFactory, 2>> const& remainders,
    std::vector<cpp_rational> const& starting_arguments,
    std::function<void(/*index=*/std::int64_t,
                       absl::StatusOr<cpp_rational>)> const& callback);

}  // namespace internal

using internal::AccurateFunction;
using internal::AccuratePolynomial;
using internal::AccuratePolynomialFactory;
using internal::ApproximateFunction;
using internal::ApproximateFunctionFactory;
using internal::GalExhaustiveMultisearch;
using internal::GalExhaustiveSearch;
using internal::StehléZimmermannSimultaneousFullSearch;
using internal::StehléZimmermannSimultaneousMultisearch;
using internal::StehléZimmermannSimultaneousSearch;
using internal::StehléZimmermannSimultaneousStreamingMultisearch;

}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia

#include "functions/accurate_table_generator_body.hpp"
