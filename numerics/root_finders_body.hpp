
#pragma once

#include "numerics/root_finders.hpp"

#include <algorithm>
#include <vector>
#include <limits>

#include "geometry/barycentre_calculator.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "numerics/double_precision.hpp"
#include "numerics/scale_b.h"

namespace principia {
namespace numerics {
namespace internal_root_finders {

using geometry::Barycentre;
using geometry::Sign;
using quantities::Abs;
using quantities::Difference;
using quantities::Sqrt;
using quantities::Square;

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
  CHECK_NE(Sign(f_lower), Sign(f_upper))
      << "\nlower: " << lower_bound << u8" ↦ " << f_lower << ", "
      << "\nupper: " << upper_bound << u8" ↦ " << f_upper;
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

// The implementation is translated from the ALGOL 60 in [Bre73], chapter 4,
// section 6, with the notation adjusted to more closely mirror the formulæ from
// the preceding sections.
template<typename Argument, typename Function>
Argument Brent(Function f,
               Argument const& lower_bound,
               Argument const& upper_bound) {
  using Value = decltype(f(lower_bound));
  Value const zero{};


  // We do not use |std::numeric_limits<double>::epsilon()|, because it is 2ϵ in
  // Brent’s notation: Brent uses ϵ = β^(1-τ) / 2 for rounded arithmetic, see
  // (2.9).
  constexpr double ϵ = ScaleB(0.5, 1 - std::numeric_limits<double>::digits);

  Argument a = lower_bound;
  Argument b = upper_bound;
  Argument c;

  Difference<Argument> d;
  Difference<Argument> e;

  Value f_a = f(a);
  Value f_b = f(b);
  Value f_c;

  if (f_a == zero) {
    return a;
  }
  if (f_b == zero) {
    return b;
  }
  CHECK_NE(Sign(f_a), Sign(f_b))
      << "\nlower: " << lower_bound << u8" ↦ " << f_a << ", "
      << "\nupper: " << upper_bound << u8" ↦ " << f_b;

  for (;;) {
    c = a;
    f_c = f_a;
    d = e = b - a;
    do {
      if (Abs(f_c) < Abs(f_b)) {
        a = b;
        b = c;
        c = a;
        f_a = f_b;
        f_b = f_c;
        f_c = f_a;
      }
      Difference<Argument> const δ = 2 * ϵ * Abs(b - Argument{});
      Difference<Argument> const m = 0.5 * (c - b);
      if (Abs(m) <= δ || f_b == zero) {
        return b;
      }
      // See if a bisection is forced.
      if (Abs(e) < δ || Abs(f_a) <= Abs(f_b)) {
        d = e = m;
      } else {
        Difference<Argument> p;
        double q;
        if (a == c) {
          // Linear interpolation.
          double const s = f_b / f_a;
          p = 2 * m * s;
          q = 1 - s;
        } else {
          // Inverse quadratic interpolation.
          double const r₁ = f_a / f_c;
          double const r₂ = f_b / f_c;
          double const r₃ = f_b / f_a;
          p = r₃ * (2 * m * r₁ * (r₁ - r₂) - (b - a) * (r₂ - 1));
          q = (r₁ - 1) * (r₂ - 1) * (r₃ - 1);
        }
        if (Sign(p).is_positive()) {
          q = -q;
        } else {
          p = -p;
        }
        if (2 * p < 3 * m * q - Abs(δ * q) && p < Abs(0.5 * e * q)) {
          e = d;
          d = p / q;
        } else {
          d = e = m;
        }
      }
      a = b;
      f_a = f_b;
      b += Abs(d) > δ ? d : δ * Sign(m);
      f_b = f(b);
    } while (Sign(f_b) != Sign(f_c));
  }
}

// See https://en.wikipedia.org/wiki/Golden-section_search for a description of
// this algorithm.
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
