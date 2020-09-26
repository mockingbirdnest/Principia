
#pragma once

#include "numerics/root_finders.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <vector>

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
using quantities::Product;
using quantities::Sqrt;
using quantities::Square;
namespace si = quantities::si;

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
                             Compare const compare) {
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

  int ev = 4;

  while (lower < lower_interior &&
         lower_interior < upper_interior &&
         upper_interior < upper) {
    Value const f_lower_min = std::min(f_lower, f_lower_interior, compare);
    Value const f_upper_min = std::min(f_upper_interior, f_upper, compare);
    if (compare(f_lower_min, f_upper_min)) {
      upper = upper_interior;
      f_upper = f_upper_interior;
      upper_interior = lower_interior;
      f_upper_interior = f_lower_interior;
      // The new lower interior point must be not be computed using the upper
      // point, lest the ratios diverge.  A very similar issue is discussed in
      // [Ove65] (who does not mention the resolution).  We use the formula from
      // the golden section step in [Bre73], chapter 5, section 8.
      lower_interior =
          upper_interior - (upper_interior - lower) * lower_interior_ratio;
      f_lower_interior = f(lower_interior);
      ++ev;
    } else {
      lower = lower_interior;
      f_lower = f_lower_interior;
      lower_interior = upper_interior;
      f_lower_interior = f_upper_interior;
      // See the remark in the other branch.
      upper_interior =
          lower_interior + (upper - lower_interior) * lower_interior_ratio;
      f_upper_interior = f(upper_interior);
      ++ev;
    }
  }
  LOG(ERROR)<<ev;

  return Barycentre<Argument, double>({lower, upper}, {1, 1});
}

// The implementation is translated from the ALGOL 60 in [Bre73], chapter 5,
// section 8.
template<typename Argument, typename Function, typename Compare>
Argument Brent(Function f,
               Argument const& lower_bound,
               Argument const& upper_bound,
               Compare compare,
               double eps) {
  using Value = decltype(f(lower_bound));

  static_assert(std::is_same_v<Compare, std::less<>> ||
                    std::is_same_v<Compare, std::greater<>>,
                "Brent’s method relies on the consistency of the order whose "
                "extremum is sought with the arithmetic operations.  For "
                "arbitrary order relations, use golden section search.");

  // The code from [Bre73] looks for a minimum; for a maximum, we look for a
  // minimum of the opposite.
  auto const minimized_f = [&f](Argument const& x) {
    if constexpr (std::is_same_v<Compare, std::greater<>>) {
      return -f(x);
    } else {
      return f(x);
    }
  };
  {
    auto const& f = minimized_f;

    // We do not use |std::numeric_limits<double>::epsilon()|, because it is 2ϵ
    // in Brent’s notation: Brent uses ϵ = β^(1-τ) / 2 for rounded arithmetic,
    // see [Bre73], chapter 4, (2.9).
    constexpr double ϵ = ScaleB(0.5, 1 - std::numeric_limits<double>::digits);
    // In order to ensure convergence, eps should be no smaller than 2ϵ, see
    // [Bre73] chapter 5, section 5.
    eps = std::max(eps, 2 * ϵ);
    // Similarly, t needs to be greater than 0, see [Bre73] chapter 5,
    // section 4.
    constexpr Difference<Argument> t =
        std::numeric_limits<double>::denorm_min() *
        si::Unit<Difference<Argument>>;

    Argument a = lower_bound;
    Argument b = upper_bound;
    constexpr double c = 2 - φ;
    Difference<Argument> d;
    Argument u;
    Argument v;
    Argument w;
    Argument x;
    Value f_u;
    Value f_v;
    Value f_w;
    Value f_x;

    v = w = x = a + c * (b - a);
    Difference<Argument> e{};
    f_v = f_w = f_x = f(x);
    int ev = 1;
    for (;;) {
      Argument const m = Barycentre<Argument, double>({a, b}, {1, 1});
      Difference<Argument> const tol = eps * Abs(x - Argument{}) + t;
      Difference<Argument> const t2 = 2 * tol;
      // Check stopping criterion.
      if (Abs(x - m) <= t2 - 0.5 * (b - a)) {
        LOG(ERROR)<<ev;
        return x;
      }
      // p = q = r = 0;
      Product<Square<Difference<Argument>>, Difference<Value>> p{};
      Product<Difference<Argument>, Difference<Value>> q{};
      if (Abs(e) > tol) {
        // Fit parabola.
        auto const r₁ = (x - w) * (f_x - f_v);
        auto const r₂ = (x - v) * (f_x - f_w);
        p = (x - v) * r₂ - (x - w) * r₁;
        q = 2 * (r₂ - r₁);
        if (Sign(q).is_positive()) {
          p = -p;
        } else {
          q = -q;
        }
      }
      // The second clause is incorrectly p < q * (a - x) in [Bre73] p.80, see
      // the errata.
      if (Abs(p) < Abs(0.5 * q * e) && p > q * (a - x) && p < q * (b - x)) {
        e = d;
        // A “parabolic interpolation” step.
        d = p / q;
        u = x + d;
        // f must not be evaluated too close to a or b.
        if (u - a < t2 || b - u < t2) {
          d = x < m ? tol : -tol;
        }
      } else {
        // A “golden section” step.
        e = (x < m ? b : a) - x;
        d = c * e;
      }
      // f must not be evaluated too close to x.
      u = x + (Abs(d) > tol ? d : tol * Sign(d));
      f_u = f(u);
      ++ev;
      // Update a, b, v, w, and x.
      if (f_u <= f_x) {
        if (u < x) {
          b = x;
        } else {
          a = x;
        }
        v = w;
        f_v = f_w;
        w = x;
        f_w = f_x;
        x = u;
        f_x = f_u;
      } else {
        if (u < x) {
          a = u;
        } else {
          b = u;
        }
        if (f_u <= f_w || w == x) {
          v = w;
          f_v = f_w;
          w = u;
          f_w = f_u;
        } else if (f_u <= f_v || v == x || v == w) {
          v = u;
          f_v = f_u;
        }
      }
    }
  }
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
