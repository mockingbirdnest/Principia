#pragma once

#include <list>
#include <type_traits>

#include "base/ranges.hpp"
#include "numerics/hermite3.hpp"

namespace principia {
namespace numerics {
namespace internal_fit_hermite_spline {

using base::Range;
using geometry::Normed;
using quantities::Derivative;

template<typename Samples,
         typename GetArgument,
         typename GetValue,
         typename GetDerivative,
         typename ErrorType>
std::list<typename Samples::const_iterator> FitHermiteSpline(
    Samples const& samples,
    GetArgument const& get_argument,
    GetValue const& get_value,
    GetDerivative const& get_derivative,
    ErrorType const& tolerance) {
  using Iterator = typename Samples::const_iterator;
  using Argument = std::remove_const_t<
      std::remove_reference_t<decltype(get_argument(*samples.begin()))>>;
  using Value = std::remove_const_t<
      std::remove_reference_t<decltype(get_value(*samples.begin()))>>;
  using Derivative1 = std::remove_const_t<
      std::remove_reference_t<decltype(get_derivative(*samples.begin()))>>;
  static_assert(
      std::is_same<Derivative1, Derivative<Value, Argument>>::value,
      "Inconsistent types for |get_argument|, |get_value|, and "
      "|get_derivative|");
  static_assert(
      std::is_same<ErrorType, typename Normed<Value>::NormType>::value,
      "|tolerance| must have the same type as a distance between values");

  std::ptrdiff_t const size = samples.end() - samples.begin();
  if (size < 3) {
    // With 0 or 1 points there is nothing to interpolate, with 2 we cannot
    // estimate the error.
    return {};
  }

  auto interpolation_error = [get_argument, get_value, get_derivative](
                                 Iterator begin, Iterator last) {
    return Hermite3<Argument, Value>(
               {get_argument(*begin), get_argument(*last)},
               {get_value(*begin), get_value(*last)},
               {get_derivative(*begin), get_derivative(*last)})
        .LInfinityError(Range(begin, last + 1), get_argument, get_value);
  };

  auto const last = samples.end() - 1;
  if (interpolation_error(samples.begin(), last) < tolerance) {
    // A single polynomial fits the entire range, so we have no way of knowing
    // whether it is the largest polynomial that will fit the range.
    return {};
  } else {
    // Look for the longest polynomial that will fit the beginning within
    // |tolerance|.
    // This will find a locally longest such polynomial, but we hope that the
    // error often grows with the number of points fitted, so that it does not
    // matter much.  In any case, the chosen polynomial does fit the beginning
    // of |samples| within |tolerance|.

    // Invariant: The Hermite interpolant on [samples.begin(), lower] is below
    // the tolerance, the Hermite interpolant on [samples.begin(), upper] is
    // above.
    auto lower = samples.begin() + 1;
    auto upper = last;
    for (;;) {
      auto const middle = lower + (upper - lower) / 2;
      if (middle == lower) {
        break;
      }
      if (interpolation_error(samples.begin(), middle) < tolerance) {
        lower = middle;
      } else {
        upper = middle;
      }
    }
    std::list<Iterator> tail = FitHermiteSpline(Range(lower, samples.end()),
                                                get_argument,
                                                get_value,
                                                get_derivative,
                                                tolerance);
    tail.push_front(lower);
    return tail;
  }
}

}  // namespace internal_fit_hermite_spline
}  // namespace numerics
}  // namespace principia
