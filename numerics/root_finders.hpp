#pragma once

#include <functional>
#include <limits>
#include <vector>

#include "absl/container/btree_set.h"
#include "base/array.hpp"
#include "numerics/scale_b.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _root_finders {
namespace internal {

using namespace principia::base::_array;
using namespace principia::numerics::_scale_b;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_named_quantities;

// Approximates a root of `f` between `lower_bound` and `upper_bound` by
// bisection.  The result is less than one ULP from a root of any continuous
// function agreeing with `f` on the values of `Argument`.
// If `f(lower_bound)` and `f(upper_bound)` are both nonzero, they must be of
// opposite signs.
template<typename Argument, typename Function>
Argument Bisect(Function f,
                Argument const& lower_bound,
                Argument const& upper_bound);

// Performs Brent’s procedure `zero` from [Bre73], chapter 4, with an absolute
// tolerance t=0.
template<typename Argument, typename Function>
Argument Brent(Function f,
               Argument const& lower_bound,
               Argument const& upper_bound);

// Alternatively applies Brent’s procedure `zero` and Brent’s procedure
// `localmin` (both from [Bre73]) to find all the zeroes over the interval by
// separating them using the extrema.  `eps` must be significantly larger than
// the default `eps` passed to `localmin`.  A value of `eps` that's too small
// will cause the algorithm to loop forever because calls to `localmin` with
// different bounds (around the same location) may return slightly different
// values for the extremum, and confuse the separation.  The default value was
// chosen experimentally: large enough to avoid that problem, and small enough
// to properly separate the zeroes.
template<typename Argument, typename Function>
absl::btree_set<Argument> DoubleBrent(
    Function f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    double eps = 1000 *
                 Sqrt(ScaleB(0.5, 1 - std::numeric_limits<double>::digits)));

// Performs a golden-section search to find a local extremum of `f` between
// `lower_bound` and `upper_bound`.
// The function searches for a minimum if compare is <, and a maximum if compare
// is >.  Arbitrary order relations are allowed; in general, this function
// searches for a value x such that compare(y, f(x)) is false for all y in some
// neighbourhood of x.
template<typename Argument, typename Function, typename Compare>
Argument GoldenSectionSearch(Function f,
                             Argument const& lower_bound,
                             Argument const& upper_bound,
                             Compare compare);

// Performs Brent’s procedure `localmin` from [Bre73], chapter 5, with an
// absolute tolerance t set to the (subnormal) smallest strictly positive value
// of `Difference<Argument>`.
// The function searches for a minimum if compare is <, and a maximum if compare
// is >.  No values of Compare other than std::less<> and std::greater<> are
// allowed.
// The default value of `eps` is √ϵ, for ϵ as defined in [Bre73], chapter 4,
// (2.9).
template<typename Argument, typename Function, typename Compare>
Argument Brent(
    Function f,
    Argument const& lower_bound,
    Argument const& upper_bound,
    Compare compare,
    double eps = Sqrt(ScaleB(0.5, 1 - std::numeric_limits<double>::digits)));

// Returns the solutions of the quadratic equation:
//   a₂ * (x - origin)² + a₁ * (x - origin) + a₀ == 0
// The result may have 0, 1 or 2 values and is sorted.
template<typename Argument, typename Value>
BoundedArray<Argument, 2> SolveQuadraticEquation(
    Argument const& origin,
    Value const& a₀,
    Derivative<Value, Argument> const& a₁,
    Derivative<Derivative<Value, Argument>, Argument> const& a₂);

}  // namespace internal

using internal::Bisect;
using internal::Brent;
using internal::DoubleBrent;
using internal::GoldenSectionSearch;
using internal::SolveQuadraticEquation;

}  // namespace _root_finders
}  // namespace numerics
}  // namespace principia

#include "numerics/root_finders_body.hpp"
