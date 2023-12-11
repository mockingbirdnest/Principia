# pragma once

#include <type_traits>

#include "numerics/чебышёв_series.hpp"

namespace principia {
namespace numerics {
namespace _approximation {
namespace internal {

using namespace principia::numerics::_чебышёв_series;

template<typename Argument, typename Function>
ЧебышёвSeries<std::invoke_result_t<Function, Argument>>
ЧебышёвPolynomialInterpolant(Function f,
                             Argument const& lower_bound,
                             Argument const& upper_bound);

}  // namespace internal

using internal::ЧебышёвPolynomialInterpolant;

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia

#include "numerics/approximation_body.hpp"
