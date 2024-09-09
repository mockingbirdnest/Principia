#pragma once

#include "numerics/fit_hermite_spline.hpp"

#include <list>
#include <type_traits>

#include "base/jthread.hpp"  // 🧙 For RETURN_IF_STOPPED.
#include "base/ranges.hpp"
#include "numerics/hermite3.hpp"

namespace principia {
namespace numerics {
namespace _fit_hermite_spline {
namespace internal {

using namespace principia::base::_ranges;
using namespace principia::numerics::_hermite3;

template<typename Value, typename Argument, typename Samples>
absl::StatusOr<std::list<typename Samples::const_iterator>> FitHermiteSpline(
    Samples const& samples,
    std::function<Argument const&(typename Samples::value_type const&)> const&
        get_argument,
    std::function<Value const&(typename Samples::value_type const&)> const&
        get_value,
    std::function<Derivative<Value, Argument> const&(
        typename Samples::value_type const&)> const& get_derivative,
    typename Hilbert<Difference<Value>>::NormType const& tolerance) {
  using Iterator = typename Samples::const_iterator;

  auto interpolation_error_is_within_tolerance =
      [get_argument, get_derivative, get_value, tolerance](
          Iterator const begin, Iterator const last) {
        return Hermite3<Value, Argument>(
                   {get_argument(*begin), get_argument(*last)},
                   {get_value(*begin), get_value(*last)},
                   {get_derivative(*begin), get_derivative(*last)})
            .LInfinityErrorIsWithin(
                Range(begin, last + 1), get_argument, get_value, tolerance);
      };

  std::list<Iterator> fit;
  if (samples.size() < 3) {
    // With 0 or 1 points there is nothing to interpolate, with 2 we cannot
    // estimate the error.
    return fit;
  }

  Iterator begin = samples.begin();
  Iterator const last = samples.end() - 1;
  while (last - begin + 1 >= 3 &&
         !interpolation_error_is_within_tolerance(begin, last)) {
    // Look for a cubic that fits the beginning within `tolerance` and
    // such the cubic fitting one more sample would not fit the samples within
    // `tolerance`.
    // Note that there may be more than one cubic satisfying this property;
    // ideally we would like to find the longest one, but this would be costly,
    // and we do not expect significant gains from this in practice.

    // Invariant: The Hermite interpolant on [begin, lower] is below the
    // tolerance, the Hermite interpolant on [begin, upper] is above.
    Iterator lower = begin + 1;
    Iterator upper = last;
    for (;;) {
      RETURN_IF_STOPPED;
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
      if (interpolation_error_is_within_tolerance(begin, middle)) {
        lower = middle;
      } else {
        upper = middle;
      }
    }
    fit.push_back(lower);

    begin = lower;
  }

  // If downsampling is not effective we'll output one iterator for each input
  // point, except at the end where we give up because we don't have enough
  // points left.
#if PRINCIPIA_MUST_ALWAYS_DOWNSAMPLE
  CHECK_LT(fit.size(), samples.size() - 2);
#endif
  return fit;
}

}  // namespace internal
}  // namespace _fit_hermite_spline
}  // namespace numerics
}  // namespace principia
