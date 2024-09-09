#pragma once

#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#include "base/not_null.hpp"
#include "numerics/polynomial_in_чебышёв_basis.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _approximation {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::numerics::_polynomial_in_чебышёв_basis;
using namespace principia::quantities::_named_quantities;

template<typename Argument, typename Function>
using Value = std::invoke_result_t<Function, Argument>;

// A function that returns true iff the interpolation interval should be split
// further.  It is only called if the `error_estimate` is larger than the
// `max_error` given to `AdaptiveЧебышёвPolynomialInterpolant`.
template<typename Value, typename Argument>
using SubdivisionPredicate = std::function<bool(
    PolynomialInЧебышёвBasis<Value, Argument> const& interpolant,
    Difference<Value> const& error_estimate)>;

// A function that returns true iff construction of the interpolants must stop.
template<typename Value, typename Argument>
using TerminationPredicate = std::function<bool(
    not_null<std::unique_ptr<PolynomialInЧебышёвBasis<Value, Argument>>>
        interpolant)>;

// Returns a Чебышёв polynomial interpolant of f over
// [lower_bound, upper_bound].  Stops if the absolute error is estimated to be
// below `max_error` or if `max_degree` has been reached.  If `error_estimate`
// is nonnull, it receives the estimate of the L∞ error.
template<int max_degree, typename Argument, typename Function>
not_null<std::unique_ptr<
    PolynomialInЧебышёвBasis<Value<Argument, Function>, Argument>>>
ЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    Difference<Value<Argument, Function>>* error_estimate = nullptr);

// Returns an ordered vector of Чебышёв polynomial interpolants of f, which
// together cover [lower_bound, upper_bound].  Subdivides the interval until the
// error is below `max_error` or `subdivide` returns false.
template<int max_degree, typename Argument, typename Function>
std::vector<not_null<std::unique_ptr<
    PolynomialInЧебышёвBasis<Value<Argument, Function>, Argument>>>>
AdaptiveЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    SubdivisionPredicate<Value<Argument, Function>, Argument> const& subdivide,
    Difference<Value<Argument, Function>>* error_estimate = nullptr);

// A streaming version of the above: as each interpolant that is below
// `max_error` is computed, it is passed to `stop`, which should return false
// if the production of interpolants should continue and true if it should stop.
template<int max_degree, typename Argument, typename Function>
void StreamingAdaptiveЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    SubdivisionPredicate<Value<Argument, Function>, Argument> const& subdivide,
    TerminationPredicate<Value<Argument, Function>, Argument> const& stop,
    Difference<Value<Argument, Function>>* error_estimate = nullptr);

}  // namespace internal

using internal::AdaptiveЧебышёвPolynomialInterpolant;
using internal::StreamingAdaptiveЧебышёвPolynomialInterpolant;
using internal::SubdivisionPredicate;
using internal::TerminationPredicate;
using internal::ЧебышёвPolynomialInterpolant;

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia

#include "numerics/approximation_body.hpp"
