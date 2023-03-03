#pragma once

#include <functional>
#include <limits>
#include <vector>

#include "base/array.hpp"
#include "numerics/scale_b.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_root_finders {

using namespace principia::base::_array;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

// Approximates a root of |f| between |lower_bound| and |upper_bound| by
// bisection.  The result is less than one ULP from a root of any continuous
// function agreeing with |f| on the values of |Argument|.
// If |f(lower_bound)| and |f(upper_bound)| are both nonzero, they must be of
// opposite signs.
template<typename Argument, typename Function>
Argument Bisect(Function f,
                Argument const& lower_bound,
                Argument const& upper_bound);

// Performs Brent’s procedure |zero| from [Bre73], chapter 4, with an absolute
// tolerance t=0.
template<typename Argument, typename Function>
Argument Brent(Function f,
               Argument const& lower_bound,
               Argument const& upper_bound);

// Performs a golden-section search to find a local extremum of |f| between
// |lower_bound| and |upper_bound|.
// The function searches for a minimum if compare is <, and a maximum if compare
// is >.  Arbitrary order relations are allowed; in general, this function
// searches for a value x such that compare(y, f(x)) is false for all y in some
// neighbourhood of x.
template<typename Argument, typename Function, typename Compare>
Argument GoldenSectionSearch(Function f,
                             Argument const& lower_bound,
                             Argument const& upper_bound,
                             Compare compare);

// Performs Brent’s procedure |localmin| from [Bre73], chapter 5, with an
// absolute tolerance t set to the (subnormal) smallest strictly positive value
// of |Difference<Argument>|.
// The function searches for a minimum if compare is <, and a maximum if compare
// is >.  No values of Compare other than std::less<> and std::greater<> are
// allowed.
// The default value of |eps| is √ϵ, for ϵ as defined in [Bre73], chapter 4,
// (2.9).
template<typename Argument, typename Function, typename Compare>
Argument Brent(
    Function f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Compare compare,
    double eps = Sqrt(ScaleB(0.5, 1 - std::numeric_limits<double>::digits)));

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
using internal_root_finders::Brent;
using internal_root_finders::GoldenSectionSearch;
using internal_root_finders::SolveQuadraticEquation;

}  // namespace numerics
}  // namespace principia

#include "numerics/root_finders_body.hpp"
