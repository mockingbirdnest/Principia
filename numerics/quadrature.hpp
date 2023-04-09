#pragma once

#include <optional>
#include <type_traits>

#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _quadrature {
namespace internal {

using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> GaussLegendre(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound);

// Computes a Clenshaw-Curtis quadrature on 2ᵖ + 1 points for successive p until
// the tolerance is satisfied.  The client controls the accuracy of the result
// using |max_relative_error| (returns when the relative error on the result is
// estimated to be less that the specified value) and |max_points| (returns when
// the number of points would exceed the specified value).
template<int initial_points = 3, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument>
AutomaticClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    std::optional<double> max_relative_error,
    std::optional<int> max_points);

// |points| must be of the form 2ᵖ + 1 for some p ∈ ℕ.  Returns the
// Clenshaw-Curtis quadrature of f with the given number of points.
template<int points, typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> ClenshawCurtis(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound);

// Computes a heuristic for the maximum number of points for an oscillating
// function.
std::optional<int> MaxPointsHeuristicsForAutomaticClenshawCurtis(
    AngularFrequency const& max_ω,
    Time const& Δt,
    int min_points_overall,
    int points_per_period);

template<typename Argument, typename Function>
Primitive<std::invoke_result_t<Function, Argument>, Argument> Midpoint(
    Function const& f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    int intervals);

}  // namespace internal

using internal::AutomaticClenshawCurtis;
using internal::GaussLegendre;
using internal::MaxPointsHeuristicsForAutomaticClenshawCurtis;
using internal::Midpoint;

}  // namespace _quadrature
}  // namespace numerics
}  // namespace principia

namespace principia::numerics {
using namespace principia::numerics::_quadrature;
}  // namespace principia::numerics

#include "numerics/quadrature_body.hpp"
