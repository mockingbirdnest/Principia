
#pragma once

#include <set>

#include "quantities/quantities.hpp"

namespace principia {

using quantities::Derivative;

namespace numerics {

// Approximates a root of |f| between |lower_bound| and |upper_bound| by
// bisection.  The result is less than one ULP from a root of any continuous
// function agreeing with |f| on the values of |Argument|.
// |f(lower_bound)| and |f(upper_bound)| must be nonzero and of opposite signs.
template<typename Argument, typename Function>
Argument Bisect(Function f,
                Argument const& lower_bound,
                Argument const& upper_bound);

// Returns the solutions of the quadratic equation:
//   a2 * (x - origin)^2 + a1 * (x - origin) + a0 == 0
// The result may have 0, 1 or 2 values.
template <typename Argument, typename Value>
std::set<Argument> SolveQuadraticEquation(
    Argument const& origin,
    Value const& a0,
    Derivative<Value, Argument> const& a1,
    Derivative<Derivative<Value, Argument>, Argument> const& a2);

}  // namespace numerics
}  // namespace principia

#include "numerics/root_finders_body.hpp"
