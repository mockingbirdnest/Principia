
#pragma once

#include "root_finders.hpp"

#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "numerics/double_precision.hpp"

namespace principia {
namespace numerics {
namespace internal_root_finders {

using geometry::Barycentre;
using geometry::Sign;
using quantities::Square;
using quantities::Sqrt;

template<typename Argument, typename Function>
Argument Bisect(Function f,
                Argument const& lower_bound,
                Argument const& upper_bound) {
  using Value = decltype(f(lower_bound));
  Value const zero{};
  Value f_upper = f(upper_bound);
  Value f_lower = f(lower_bound);
  if (f_upper == zero) {
    return upper_bound;
  }
  if (f_lower == zero) {
    return lower_bound;
  }
  CHECK(f_lower > zero && zero > f_upper || f_lower < zero && zero < f_upper)
      << "\nlower: " << lower_bound << " :-> " << f_lower << ", "
      << "\nupper: " << upper_bound << " :-> " << f_upper;
  Argument lower = lower_bound;
  Argument upper = upper_bound;
  for (;;) {
    Argument const middle =
        Barycentre<Argument, double>({lower, upper}, {1, 1});
    // The size of the interval has reached one ULP.
    if (middle == lower || middle == upper) {
      return middle;
    }
    Value const f_middle = f(middle);
    if (f_middle == zero) {
      return middle;
    } else if (Sign(f_middle) == Sign(f_upper)) {
      upper = middle;
    } else {
      lower = middle;
    }
  }
}

//https://en.wikipedia.org/wiki/Golden-section_search
template<typename Argument, typename Function, typename Compare>
Argument GoldenSectionSearch(Function f,
                             Argument const& lower_bound,
                             Argument const& upper_bound,
                             Compare const comp) {
  static constexpr double lower_interior_ratio = 2 - φ;
  static constexpr double upper_interior_ratio = φ - 1;
  using Value = decltype(f(lower_bound));

  Argument upper = upper_bound;
  Value f_upper = f(upper);

  Argument lower = lower_bound;
  Value f_lower = f(lower);

  Argument lower_interior = Barycentre<Argument, double>(
      {lower, upper}, {upper_interior_ratio, lower_interior_ratio});
  Value f_lower_interior = f(lower_interior);

  Argument upper_interior = Barycentre<Argument, double>(
      {lower, upper}, {lower_interior_ratio, upper_interior_ratio});
  Value f_upper_interior = f(upper_interior);

  while (lower < lower_interior &&
         lower_interior < upper_interior &&
         upper_interior < upper) {
    Value const f_lower_min = std::min(f_lower, f_lower_interior, comp);
    Value const f_upper_min = std::min(f_upper_interior, f_upper, comp);
    if (comp(f_lower_min, f_upper_min)) {
      upper = upper_interior;
      f_upper = f_upper_interior;
      upper_interior = lower_interior;
      f_upper_interior = f_lower_interior;
      lower_interior = Barycentre<Argument, double>(
          {lower, upper}, {upper_interior_ratio, lower_interior_ratio});
      f_lower_interior = f(lower_interior);
    } else {
      lower = lower_interior;
      f_lower = f_lower_interior;
      lower_interior = upper_interior;
      f_lower_interior = f_upper_interior;
      upper_interior = Barycentre<Argument, double>(
          {lower, upper}, {lower_interior_ratio, upper_interior_ratio});
      f_upper_interior = f(upper_interior);
    }
  }

  return Barycentre<Argument, double>({lower, upper}, {1, 1});
}

template<typename Argument, typename Value>
BoundedArray<Argument, 2> SolveQuadraticEquation(
    Argument const& origin,
    Value const& a0,
    Derivative<Value, Argument> const& a1,
    Derivative<Derivative<Value, Argument>, Argument> const& a2) {
  using Derivative1 = Derivative<Value, Argument>;
  using Discriminant = Square<Derivative1>;

  // This algorithm is after section 1.8 of Accuracy and Stability of Numerical
  // Algorithms, Second Edition, Higham, ISBN 0-89871-521-0.

  static Discriminant const discriminant_zero{};

  // Use double precision for the discriminant because there can be
  // cancellations.  Higham mentions that it is necessary "to use extended
  // precision (or some trick tantamount to the use of extended precision)."
  DoublePrecision<Discriminant> discriminant = TwoSum(a1 * a1, -4.0 * a0 * a2);

  if (discriminant.value == discriminant_zero &&
      discriminant.error == discriminant_zero) {
    // One solution.
    return {origin - 0.5 * a1 / a2};
  } else if (discriminant.value < discriminant_zero ||
             (discriminant.value == discriminant_zero &&
              discriminant.error < discriminant_zero)) {
    // No solution.
    return {};
  } else {
    // Two solutions.  Compute the numerator of the larger one.
    Derivative1 numerator;
    static Derivative1 derivative_zero{};
    if (a1 > derivative_zero) {
      numerator = -a1 - Sqrt(discriminant.value + discriminant.error);
    } else {
      numerator = -a1 + Sqrt(discriminant.value + discriminant.error);
    }
    auto const solution1 = origin + numerator / (2.0 * a2);
    auto const solution2 = origin + (2.0 * a0) / numerator;
    if (solution1 < solution2) {
      return {solution1, solution2};
    } else {
      return {solution2, solution1};
    }
  }
}

}  // namespace internal_root_finders
}  // namespace numerics
}  // namespace principia
