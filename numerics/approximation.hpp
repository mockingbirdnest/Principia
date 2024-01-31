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

// A function that returns true iff the interpolation interval should be split
// further.  It is only called if the |error_estimate| is larger than the
// |max_error| given to |AdaptiveЧебышёвPolynomialInterpolant|.
template<typename Value, typename Argument>
using SubdivisionPredicate =
    std::function<bool(ЧебышёвSeries<Value, Argument> const& interpolant,
                       Difference<Value> const& error_estimate)>;

// A function that returns true if construction of the interpolants
// must proceed.
template<typename Value, typename Argument>
using ProceedPredicate =
    std::function<bool(ЧебышёвSeries<Value, Argument> interpolant)>;

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

// Returns an ordered vector of Чебышёв polynomial interpolants of f, which
// together cover [lower_bound, upper_bound].  Subdivides the interval until the
// error is below |max_error| or |subdivide| returns false.
template<int max_degree, typename Argument, typename Function>
std::vector<ЧебышёвSeries<Value<Argument, Function>, Argument>>
AdaptiveЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    SubdivisionPredicate<Value<Argument, Function>, Argument> const& subdivide,
    Difference<Value<Argument, Function>>* error_estimate = nullptr);

// A streaming version of the above: as each interpolant that is below
// |max_error| is computed, it is passed to |proceed|, which should return false
// if the production of interpolants should continue and false if it should
// stop.
template<int max_degree, typename Argument, typename Function>
void StreamingAdaptiveЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Argument, Function>> const& max_error,
    SubdivisionPredicate<Value<Argument, Function>, Argument> const& subdivide,
    ProceedPredicate<Value<Argument, Function>, Argument> const& proceed,
    Difference<Value<Argument, Function>>* error_estimate = nullptr);

}  // namespace internal

using internal::AdaptiveЧебышёвPolynomialInterpolant;
using internal::StreamingAdaptiveЧебышёвPolynomialInterpolant;
using internal::SubdivisionPredicate;
using internal::ЧебышёвPolynomialInterpolant;

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia

#include "numerics/approximation_body.hpp"
