#pragma once

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

// Returns a Чебышёв polynomial interpolant of f over
// [lower_bound, upper_bound].  Stops if the absolute error is estimated to be
// below |max_error| or if |max_degree| has been reached.  If |error_estimate|
// is nonnull, it receives the estimate of the L∞ error.
template<int max_degree, typename Argument, typename Function>
ЧебышёвSeries<Value<Argument, Function>, Argument> ЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    Difference<Value<Argument, Function>>* error_estimate = nullptr);

// Returns an ordered vector of Чебышёв polynomial interpolants of f over, which
// together cover [lower_bound, upper_bound].  Subdivides the interval so that
// the degree of each approximant doesn't exceed |max_degree|.  The final
// (estimated) absolute error is guaranteed to be below |max_error|.
template<int max_degree, typename Argument, typename Function>
std::vector<ЧебышёвSeries<Value<Argument, Function>, Argument>>
AdaptiveЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    Difference<Value<Argument, Function>>* error_estimate = nullptr);

}  // namespace internal

using internal::AdaptiveЧебышёвPolynomialInterpolant;
using internal::ЧебышёвPolynomialInterpolant;

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia

#include "numerics/approximation_body.hpp"
