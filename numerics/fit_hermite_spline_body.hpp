
#pragma once

#include <list>
#include <type_traits>

#include "base/ranges.hpp"
#include "numerics/hermite3.hpp"

namespace principia {
namespace numerics {
namespace internal_fit_hermite_spline {

using base::Range;

template<typename Argument, typename Value, typename Samples>
std::list<typename Samples::const_iterator> FitHermiteSpline(
    Samples const& samples,
    std::function<Argument const&(typename Samples::value_type const&)> const&
        get_argument,
    std::function<Value const&(typename Samples::value_type const&)> const&
        get_value,
    std::function<Derivative<Value, Argument> const&(
        typename Samples::value_type const&)> const& get_derivative,
    typename Hilbert<Difference<Value>>::NormType const& tolerance) {
  using Iterator = typename Samples::const_iterator;

  auto interpolation_error = [get_argument, get_derivative, get_value](
                                 Iterator begin, Iterator last) {
    return Hermite3<Argument, Value>(
               {get_argument(*begin), get_argument(*last)},
               {get_value(*begin), get_value(*last)},
               {get_derivative(*begin), get_derivative(*last)})
        .LInfinityError(Range(begin, last + 1), get_argument, get_value);
  };

  std::list<Iterator> tail;
  if (samples.size() < 3) {
    // With 0 or 1 points there is nothing to interpolate, with 2 we cannot
    // estimate the error.
    return tail;
  }

  Iterator begin = samples.begin();
  Iterator const last = samples.end() - 1;
  while (last - begin + 1 >= 3 &&
         interpolation_error(begin, last) >= tolerance) {
    // Look for a cubic that fits the beginning within |tolerance| and
    // such the cubic fitting one more sample would not fit the samples within
    // |tolerance|.
    // Note that there may be more than one cubic satisfying this property;
    // ideally we would like to find the longest one, but this would be costly,
    // and we do not expect significant gains from this in practice.

    // Invariant: The Hermite interpolant on [begin, lower] is below the
    // tolerance, the Hermite interpolant on [begin, upper] is above.
    Iterator lower = begin + 1;
    Iterator upper = last;
    for (;;) {
      auto const middle = lower + (upper - lower) / 2;
      // Note that lower ≤ middle ≤ upper.
      // If middle - lower > 0, upper - lower > 0,
      // therefore (upper - lower) / 2  < upper - lower, thus
      // middle < upper.  It follows that upper - lower strictly decreases in
      // this iteration, since we assign middle to either lower or upper.
      // We exit when middle == lower, so the algorithm terminates.
      if (middle == lower) {
        break;
      }
      if (interpolation_error(begin, middle) < tolerance) {
        lower = middle;
      } else {
        upper = middle;
      }
    }
    tail.push_back(lower);

    begin = lower;
  }

  // If downsampling is not effective we'll output one iterator for each input
  // point, except at the end where we give up because we don't have enough
  // points left.
  CHECK_LT(tail.size(), samples.size() - 2);
  return tail;
}

}  // namespace internal_fit_hermite_spline
}  // namespace numerics
}  // namespace principia
