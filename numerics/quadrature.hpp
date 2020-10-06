
#pragma once

#include <type_traits>

#include "geometry/hilbert.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace quadrature {
namespace internal_quadrature {

using geometry::Hilbert;
using quantities::Primitive;

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> GaussLegendre(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound);

template<typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
AutomaticClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    typename Hilbert<Primitive<std::invoke_result_t<Function, Argument>,
                               Argument>>::NormType absolute_tolerance = {},
    double relative_tolerance = 0x1p-20);

template<typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> ClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::int64_t points);

template<typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> Midpoint(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    int intervals);

}  // namespace internal_quadrature

using internal_quadrature::AutomaticClenshawCurtis;
using internal_quadrature::GaussLegendre;
using internal_quadrature::Midpoint;

}  // namespace quadrature
}  // namespace numerics
}  // namespace principia

#include "numerics/quadrature_body.hpp"
