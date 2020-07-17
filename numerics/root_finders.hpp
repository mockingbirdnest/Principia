
#pragma once

#include <vector>

#include "base/array.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_root_finders {

using base::BoundedArray;
using quantities::Derivative;

// Approximates a root of |f| between |lower_bound| and |upper_bound| by
// bisection.  The result is less than one ULP from a root of any continuous
// function agreeing with |f| on the values of |Argument|.
// If |f(lower_bound)| and |f(upper_bound)| are both nonzero, they must be of
// opposite signs.
template<typename Argument, typename Function>
Argument Bisect(Function f,
                Argument const& lower_bound,
                Argument const& upper_bound);

//TODO(phl):comment
template<typename Argument, typename Compare, typename Function>
Argument GoldenSectionSearch(Function f,
                             Argument const& lower_bound,
                             Argument const& upper_bound);

// Returns the solutions of the quadratic equation:
//   a2 * (x - origin)^2 + a1 * (x - origin) + a0 == 0
// The result may have 0, 1 or 2 values and is sorted.
template<typename Argument, typename Value>
BoundedArray<Argument, 2> SolveQuadraticEquation(
    Argument const& origin,
    Value const& a0,
    Derivative<Value, Argument> const& a1,
    Derivative<Derivative<Value, Argument>, Argument> const& a2);

}  // namespace internal_root_finders

using internal_root_finders::Bisect;
using internal_root_finders::SolveQuadraticEquation;

}  // namespace numerics
}  // namespace principia

#include "numerics/root_finders_body.hpp"
