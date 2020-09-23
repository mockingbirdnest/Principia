
#include "numerics/root_finders.hpp"

#include <functional>
#include <set>
#include <vector>
#include <limits>

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using geometry::Instant;
using quantities::Acceleration;
using quantities::Length;
using quantities::Pow;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::IsEmpty;
using ::testing::Le;
namespace si = quantities::si;

namespace numerics {

class RootFindersTest : public ::testing::Test {};

// Solving Δt * Δt == n.
TEST_F(RootFindersTest, SquareRoots) {
  Instant const t_0;
  Instant const t_max = t_0 + 10 * Second;
  Length const n_max = Pow<2>(t_max - t_0) * si::Unit<Acceleration>;
  for (Length n = 1 * Metre; n < n_max; n += 1 * Metre) {
    int evaluations = 0;
    auto const equation = [t_0, n, &evaluations](Instant const& t) {
      ++evaluations;
      return Pow<2>(t - t_0) * si::Unit<Acceleration> - n;
    };
    EXPECT_THAT(Bisect(equation, t_0, t_max) - t_0,
                AlmostEquals(Sqrt(n / si::Unit<Acceleration>), 0, 1));
    if (n == 25 * Metre) {
      EXPECT_EQ(3, evaluations);
    } else {
      EXPECT_THAT(evaluations, AllOf(Ge(49), Le(58)));
    }

    evaluations = 0;
    EXPECT_THAT(Brent(equation, t_0, t_max) - t_0,
                AlmostEquals(Sqrt(n / si::Unit<Acceleration>), 0, 1));
    EXPECT_THAT(evaluations, AllOf(Ge(7), Le(15)));
  }
}

TEST_F(RootFindersTest, WilkinsGuFunction) {
  // This test constructs a pathological function for which Brent is quadratic.
  // The construction follows [WG13] and the notation therein.
  // We will search for a root of the constructed function over [a, b].
  double const a = 1;
  double const b = 2;

  // We first construct the sequence of points X at which Brent will be made to
  // evaluate the function; see section 4.1.
  std::vector<double> X;
  // The indices i such that Brent will perform a bisection at X[i].
  std::set<int> expected_bisections;

  // This tolerance is not the same as the one in our implementation of Brent:
  // this is because that one is overly close to the machine precision, so that
  // attempting to perform as many interpolation steps as theoretically possible
  // before the bisection leads to early bisection in practice, because the
  // interval shrinkage becomes too great as the convergent values get
  // discretized.
  double const δ = 0x1p-50;
  // The factor p theoretically only needs to be larger than √2, but, for
  // similar reasons, bringing it closer (e.g., 17/12) requires increasing δ,
  // which overall reduces the number of evaluations required.
  double const p = 1.5;

  X.push_back(b);
  int bisections = 0;
  for (; a < X.back(); ++bisections) {
    // k is the number of interpolation steps on the interval [a, X.back()]
    // before a bisection happens again.
    int const k = std::round(std::log2((X.back() - a) / δ));
    for (int j = 1; j <= k; ++j) {
      X.push_back(X.back() - std::pow(p, k - j) * δ);
    }
    expected_bisections.insert(X.size() - 1);
    X.push_back(a + (X.back() - a) / 2);
  }

  // [WG13], equation (5).
  auto const s = [](double xa, double xc, double xb, double fa) {
    return fa * (xb - xc) / (xa - xc);
  };

  // [WG13], equations (7) and (9).
  auto const q =
      [](double xa, double xd, double xc, double xb, double fa, double fb) {
        double const α = (xa - xd) * fb + (xd - xb) * fa;
        double const β = (xb - xd) * Pow<2>(fa) + (xd - xa) * Pow<2>(fb);
        double const γ = fa * fb * (fb - fa) * (xc - xd);
        return -2 * γ / (β + Sqrt(Pow<2>(β) - 4 * α * γ));
      };

  // Capture everything but the function itself by copy so we have some
  // confidence that the function is pure.
  // We capture f because it is defined recursively.
  // If non-null, |*evaluations| is incremented when the function is called.
  // If |expect_brent_calls| is true, this function checks that x is an element
  // of the expected sequence X.
  std::function<double(double, int*, bool)> f =
      [&f, s, q, a, b, X, expected_bisections](
          double const x,
          int* const evaluations,
          bool const expect_brent_calls) -> double {
    double const f_a = -100;
    if (evaluations != nullptr) {
      ++*evaluations;
    }
    if (x == a) {
      return f_a;
    }
    if (x == b) {
      return s(a, X[1], b, f_a);
    }
    // [WG13] define f only on X. We extend it as a step function, returning the
    // value for the element of X nearest to the given x.
    int j;
    double min_Δx = std::numeric_limits<double>::infinity();
    for (int i = 1; i < X.size() - 1; ++i) {
      double const Δx = std::abs(x - X[i]);
      if (Δx < min_Δx) {
        min_Δx = Δx;
        j = i;
      }
    }
    if (expect_brent_calls) {
      EXPECT_THAT(x, AlmostEquals(X[j], 0)) << j;
    }

    if (expected_bisections.count(j) != 0) {
      return 100;
    }
    // Note that the recursion does not count as an evaluation, and is not
    // subject to expectations.
    return q(a,
             X[j + 1],
             X[j],
             X[j - 1],
             f_a,
             f(X[j - 1], nullptr, /*expect_brent_calls=*/false));
  };

  int evaluations = 0;
  EXPECT_THAT(Brent(
                  [&f, &evaluations](double x) {
                    return f(x, &evaluations, /*expect_brent_calls=*/true);
                  },
                  a,
                  b),
      AlmostEquals(a, 2));
  EXPECT_THAT(evaluations, Eq(1304));

  evaluations = 0;
  EXPECT_THAT(Bisect(
                  [&f, &evaluations](double x) {
                    return f(x, &evaluations, /*expect_brent_calls=*/false);
                  },
                  a,
                  b),
              AlmostEquals(a, 0));
  EXPECT_THAT(evaluations, Eq(54));
}

TEST_F(RootFindersTest, GoldenSectionSearch) {
  Instant const t_0;
  auto sin = [t_0](Instant const& t) {
    return Sin((t - t_0) * Radian / Second);
  };

  // Minimum.
  for (int l = 16; l <= 47; ++l) {
    for (int u = 48; u <= 62; ++u) {
      // The result is not overly precise because near its maximum a the
      // function is:
      //   f(a) + fʺ(a) (x - a) / 2 + o((x - a)²)
      // The second order term vanishes when x and a match on the first 26
      // leading bits (roughly).
      EXPECT_THAT(
          GoldenSectionSearch(sin,
                              t_0 + l * 0.1 * Second,
                              t_0 + u * 0.1 * Second),
          AlmostEquals(t_0 + 3 * π / 2 * Second, 11'863'280, 11'863'284));
    }
  }

  // Maximum.
  EXPECT_THAT(
      (GoldenSectionSearch<Instant, decltype(sin), std::greater<double>>)(
          sin,
          t_0 + 1.5 * Second,
          t_0 + 1.6 * Second),
      AlmostEquals(t_0 + π / 2 * Second, 47453132));

  // A big interval will yield a semi-random minimum.
  EXPECT_THAT(GoldenSectionSearch(sin,
                                  t_0 - 100 * Second,
                                  t_0 + 666 * Second),
              AlmostEquals(t_0 + 119 * π / 2 * Second, 370'728));
}

TEST_F(RootFindersTest, QuadraticEquations) {
  // Golden ratio.
  auto const s1 = SolveQuadraticEquation(0.0, -1.0, -1.0, 1.0);
  EXPECT_THAT(s1,
              ElementsAre(AlmostEquals((1 - sqrt(5)) / 2, 1),
                          AlmostEquals((1 + sqrt(5)) / 2, 0)));

  // No solutions.
  auto const s2 = SolveQuadraticEquation(0.0, 1.0, 0.0, 1.0);
  EXPECT_THAT(s2, IsEmpty());

  // One solution.
  auto const s3 = SolveQuadraticEquation(0.0, 1.0, 2.0, 1.0);
  EXPECT_THAT(s3, ElementsAre(-1.0));

  // An ill-conditioned system.  I fart in its general direction.  If done
  // naively, this yields {-100032., -99968.4} according to Mathematica.
  auto const s4 = SolveQuadraticEquation(0.0,
                                         1.0000001e25,
                                         2.0000003e20,
                                         1.0000001e15);
  EXPECT_THAT(s4,
              ElementsAre(AlmostEquals(-100031.62777541532972762902, 66),
                          AlmostEquals(-99968.38222458367027247098, 65)));

  // A typed system.
  Instant const t0;
  auto const s5 = SolveQuadraticEquation(
      t0, 1.0 * Metre, 2.0 * Metre / Second, 1.0 * Metre / Second / Second);
  EXPECT_THAT(s5, ElementsAre(t0 - 1.0 * Second));
}

}  // namespace numerics
}  // namespace principia
