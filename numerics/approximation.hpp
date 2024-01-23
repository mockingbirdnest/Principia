#pragma once

#include <functional>
#include <type_traits>
#include <vector>

#include "numerics/чебышёв_series.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _approximation {
namespace internal {

using namespace principia::numerics::_чебышёв_series;
using namespace principia::quantities::_named_quantities;

template<typename Argument, typename Function>
using Value = std::invoke_result_t<Function, Argument>;

// A function that returns true iff the interpolant search should terminate
// because it has reached an acceptable solution.  A simple predicate can just
// check if |error_estimate| is below some threshold.
template<typename Value, typename Argument>
using TerminationPredicate =
    std::function<bool(ЧебышёвSeries<Value, Argument> const& interpolant,
                       Difference<Value> const& error_estimate)>;

// Returns a Чебышёв polynomial interpolant of f over
// [lower_bound, upper_bound].  Stops if the termination predicate |done|
// returns true or if |max_degree| has been reached.  If |error_estimate| is
// nonnull, it receives the estimate of the L∞ error.
template<int max_degree, typename Argument, typename Function>
ЧебышёвSeries<Value<Argument, Function>, Argument> ЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    TerminationPredicate<Value<Argument, Function>, Argument> const& done,
    Difference<Value<Argument, Function>>* error_estimate = nullptr);

// Returns an ordered vector of Чебышёв polynomial interpolants of f, which
// together cover [lower_bound, upper_bound].  Subdivides the interval so that
// the degree of each approximant doesn't exceed |max_degree|.
template<int max_degree, typename Argument, typename Function>
std::vector<ЧебышёвSeries<Value<Argument, Function>, Argument>>
AdaptiveЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    TerminationPredicate<Value<Argument, Function>, Argument> const& done,
    Difference<Value<Argument, Function>>* error_estimate = nullptr);

}  // namespace internal

using internal::AdaptiveЧебышёвPolynomialInterpolant;
using internal::ЧебышёвPolynomialInterpolant;
using internal::TerminationPredicate;///REMOVE

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia

#include "numerics/approximation_body.hpp"
