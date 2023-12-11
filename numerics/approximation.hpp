#pragma once

#include <type_traits>

#include "numerics/чебышёв_series.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _approximation {
namespace internal {

using namespace principia::numerics::_чебышёв_series;
using namespace principia::quantities::_named_quantities;

template<typename Function, typename Argument>
using Value = std::invoke_result_t<Function, Argument>;

template<typename Argument, typename Function>
ЧебышёвSeries<Value<Function, Argument>> ЧебышёвPolynomialInterpolant(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Difference<Value<Function, Argument>> const& max_error);

}  // namespace internal

using internal::ЧебышёвPolynomialInterpolant;

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia

#include "numerics/approximation_body.hpp"
