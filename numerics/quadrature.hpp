
#pragma once

#include "quantities/named_quantities.hpp"

#include <type_traits>

namespace principia {
namespace numerics {
namespace quadrature {
namespace internal_quadrature {

using quantities::Primitive;

template<typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> Midpoint(
    Function const& function,
    Argument const& lower_bound,
    Argument const& upper_bound,
    int intervals);

}  // namespace internal_quadrature

using internal_quadrature::Midpoint;

}  // namespace quadrature
}  // namespace numerics
}  // namespace principia

#include "numerics/quadrature_body.hpp"