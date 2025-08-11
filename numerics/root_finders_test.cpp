#include "numerics/root_finders.hpp"

#include <functional>
#include <set>
#include <vector>
#include <limits>

#include "absl/base/casts.h"
#include "geometry/instant.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "numerics/scale_b.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::IsEmpty;
using ::testing::Le;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_point;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::numerics::_root_finders;
using namespace principia::numerics::_scale_b;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

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
  // evaluate the function; see section 4.1. We denote the sequence of points by
  // X as in section 4.2; however, we number them upward starting at 0, instead
  // of downward starting at m = |X|.
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
  while (a < X.back()) {
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
  auto const s =
      [](double const xa, double const xc, double const xb, double const fa) {
        return fa * (xb - xc) / (xa - xc);
      };

  // [WG13], equations (7) and (9).
  auto const q = [](double const xa,
                    double const xd,
                    double const xc,
                    double const xb,
                    double const fa,
                    double const fb) {
    double const α = (xa - xd) * fb + (xd - xb) * fa;
    double const β = (xb - xd) * Pow<2>(fa) + (xd - xa) * Pow<2>(fb);
    double const γ = fa * fb * (fb - fa) * (xc - xd);
    return -2 * γ / (β + Sqrt(Pow<2>(β) - 4 * α * γ));
  };

  // Capture everything but the function itself by copy so we have some
  // confidence that the function is pure.
  // We capture f because it is defined recursively.
  // If non-null, `*evaluations` is incremented when the function is called.
  // If `expect_brent_calls` is true, this function checks that x is an element
  // of the expected sequence X.
  // For x in X, this function returns f(x) as defined at the end of
  // section 4.2.
  std::function<double(double, int*, bool)> f =
      [&f, s, q, a, b, X, expected_bisections](
          double const x,
          int* const evaluations,
          bool const expect_brent_calls) -> double {
    // “For the initial interval to be valid, we choose f (a) = ϵ for some
    // ϵ < 0.”
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
    int k = 0;
    double min_Δx = std::numeric_limits<double>::infinity();
    for (int i = 1; i < X.size() - 1; ++i) {
      double const Δx = std::abs(x - X[i]);
      if (Δx < min_Δx) {
        min_Δx = Δx;
        k = i;
      }
    }
    CHECK_LE(1, k);
    if (expect_brent_calls) {
      EXPECT_THAT(x, AlmostEquals(X[k], 0)) << k;
    }

    if (expected_bisections.count(k) != 0) {
      // “Since we are bisecting the function regardless of its value at X[k] ,
      // we may choose an arbitrary positive value for f (X[k]), such as 100.”
      return 100;
    }
    // Note that the recursion does not count as an evaluation, and is not
    // subject to expectations.
    return q(a, X[k + 1], X[k], X[k - 1], f_a, f(X[k - 1],
                                                 /*evaluations=*/nullptr,
                                                 /*expect_brent_calls=*/false));
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

TEST_F(RootFindersTest, RootNear0) {
  int evaluations;
  auto const f = [&evaluations](double const x) {
    ++evaluations;
    return Sin(x * Radian) - std::numeric_limits<double>::denorm_min();
  };
  evaluations = 0;
  EXPECT_THAT(Bisect(f, -1.0, 2.0),
              AlmostEquals(std::numeric_limits<double>::denorm_min(), 0));
  EXPECT_THAT(evaluations, Eq(1078));
  evaluations = 0;
  EXPECT_THAT(Brent(f, -1.0, 2.0),
              AlmostEquals(std::numeric_limits<double>::denorm_min(), 0));
  EXPECT_THAT(evaluations, Eq(10));
}

TEST_F(RootFindersTest, ManyRoots) {
  Instant const t0;
  int evaluations;
  auto const f = [&evaluations, t0](Instant const t) {
    ++evaluations;
    return Sin((t - t0) * Radian / Second);
  };
  evaluations = 0;
  EXPECT_THAT(DoubleBrent(f, t0 + 1.0 * Second, t0 + 30.0 * Second),
              ElementsAre(AlmostEquals(t0 + π * Second, 0),
                          AlmostEquals(t0 + 2 * π * Second, 0),
                          AlmostEquals(t0 + 3 * π * Second, 0),
                          AlmostEquals(t0 + 4 * π * Second, 0),
                          AlmostEquals(t0 + 5 * π * Second, 1),
                          AlmostEquals(t0 + 6 * π * Second, 0),
                          AlmostEquals(t0 + 7 * π * Second, 0),
                          AlmostEquals(t0 + 8 * π * Second, 0),
                          AlmostEquals(t0 + 9 * π * Second, 0)));
  EXPECT_THAT(evaluations, Eq(1180));
  evaluations = 0;
  EXPECT_THAT(DoubleBrent(f, t0 + π * Second, t0 + 10.0 * Second),
              ElementsAre(AlmostEquals(t0 + π * Second, 0),
                          AlmostEquals(t0 + 2 * π * Second, 0),
                          AlmostEquals(t0 + 3 * π * Second, 0)));
  EXPECT_THAT(evaluations, Eq(332));
  evaluations = 0;
  EXPECT_THAT(DoubleBrent(f,
                          t0 + π * Second,
                          t0 + std::nextafter(2 * π, 10) * Second),
              ElementsAre(AlmostEquals(t0 + π * Second, 0),
                          AlmostEquals(t0 + 2 * π * Second, 1)));
  EXPECT_THAT(evaluations, Eq(158));
  evaluations = 0;
  EXPECT_THAT(DoubleBrent(f, t0 + π * Second, t0 + 3 * π / 2 * Second),
              ElementsAre(AlmostEquals(t0 + π * Second, 0)));
  EXPECT_THAT(evaluations, Eq(76));
}

TEST_F(RootFindersTest, RootAt0) {
  int evaluations;
  auto id = [&evaluations](double const x) {
    ++evaluations;
    return x;
  };

  evaluations = 0;
  EXPECT_THAT(Bisect(id, -1.0, 2.0), AlmostEquals(0, 0));
  EXPECT_THAT(evaluations, Eq(1077));

  evaluations = 0;
  EXPECT_THAT(Brent(id, -1.0, 2.0), AlmostEquals(0, 0));
  EXPECT_THAT(evaluations, Eq(3));
}

TEST_F(RootFindersTest, MinimumNear0) {
  int evaluations;
  auto const f = [&evaluations](double const x) {
    ++evaluations;
    return Abs(x - std::numeric_limits<double>::denorm_min());
  };
  evaluations = 0;
  EXPECT_THAT(GoldenSectionSearch(f, -1.0, 1.0, std::less<>()),
              AlmostEquals(std::numeric_limits<double>::denorm_min(), 0));
  EXPECT_THAT(evaluations, Eq(1551));
  // This fails to terminate if t=0.
  evaluations = 0;
  EXPECT_THAT(Brent(f, -1.0, 1.0, std::less<>(), /*eps=*/0),
              AlmostEquals(std::numeric_limits<double>::denorm_min(), 0));
  EXPECT_THAT(evaluations, Eq(1327));
}

TEST_F(RootFindersTest, MinimumAt0) {
  int evaluations;
  auto abs = [&evaluations](double const x) {
    ++evaluations;
    return Abs(x);
  };

  // This search takes a lot of steps, so that it fails if the computation of
  // the new interior points is unstable, e.g., if they are computed as a
  // barycentre of the lower and upper points.
  evaluations = 0;
  EXPECT_THAT(GoldenSectionSearch(abs, -1.0, 1.0, std::less<>()),
              AlmostEquals(0, 0));
  EXPECT_THAT(evaluations, Eq(1551));

  evaluations = 0;
  EXPECT_THAT(Brent(abs, -1.0, 1.0, std::less<>(), 0), AlmostEquals(0, 1));
  EXPECT_THAT(evaluations, Eq(1327));
}

TEST_F(RootFindersTest, SharpMinimum) {
  int evaluations = 0;
  constexpr Instant t0;
  constexpr Point<Entropy> s0;
  // Whereas root finding requires the result type to have a 0, minimization
  // works on oriented one-dimensional affine spaces.
  auto f = [s0, t0, &evaluations](Instant const x) -> Point<Entropy> {
    ++evaluations;
    return s0 + Abs(x - (t0 + 1 * Second)) * (1 * Watt / Kelvin);
  };

  evaluations = 0;
  EXPECT_THAT(
      GoldenSectionSearch(f, t0 - π * Second, t0 + π * Second, std::less<>()),
      AlmostEquals(t0 + 1 * Second, 0));
  EXPECT_THAT(evaluations, Eq(82));

  evaluations = 0;
  EXPECT_THAT(
      Brent(f, t0 - π * Second, t0 + π * Second, std::less<>(), /*eps=*/0),
      AlmostEquals(t0 + 1 * Second, 0));
  EXPECT_THAT(evaluations, Eq(51));
}

TEST_F(RootFindersTest, SmoothMaximum) {
  int evaluations;
  // The composition of the 16th degree Taylor series for the cosine with the
  // polynomial x ↦ 3(x-1); this function approximates cos(3(x-1)) near 1, where
  // it has a local maximum.
  PolynomialInMonomialBasis<double, double, 16> const
      p({-4059064033.0 / 4100096000,
         759417921.0 / 1793792000,
         3196519569.0 / 717516800,
         -4649859.0 / 7321600,
         -7526709.0 / 2252800,
         5623263.0 / 19712000,
         3596319.0 / 3584000,
         -22599.0 / 358400,
         -3183543.0 / 20070400,
         12393.0 / 2508800,
         9477.0 / 512000,
         -6561.0 / 2816000,
         -2187.0 / 15769600,
         -19683.0 / 51251200,
         19683.0 / 102502400,
         -59049.0 / 1793792000,
         59049.0 / 28700672000},
        with_evaluator<EstrinWithoutFMA>);
  auto const f = [&evaluations, &p](double const x) {
    ++evaluations;
    return p(x);
  };

  evaluations = 0;
  EXPECT_THAT(GoldenSectionSearch(f, 0.0, π / 2, std::greater<>()),
              AlmostEquals(1, 22'535'664));
  EXPECT_THAT(evaluations, Eq(79));
  evaluations = 0;
  EXPECT_THAT(GoldenSectionSearch(f, π / 7, 9 * π / 14, std::greater<>()),
              AlmostEquals(1, 9'447'534));
  EXPECT_THAT(evaluations, Eq(80));

  constexpr double ϵ = ScaleB(0.5, 1 - std::numeric_limits<double>::digits);
  constexpr double ϵ² = ϵ * ϵ;
  constexpr double ϵ³ = ϵ² * ϵ;
  constexpr double ϵ⁵ = ϵ³ * ϵ²;

  // Locate a maximum of the computed function with full precision, starting
  // from two different intervals.
  double eps = 2 * ϵ;
  evaluations = 0;
  EXPECT_THAT(Brent(f, 0.0, π / 2, std::greater<>(), eps),
              AlmostEquals(1, 39'406'981));
  EXPECT_THAT(evaluations, Eq(37));
  evaluations = 0;
  EXPECT_THAT(Brent(f, π / 7, 9 * π / 14, std::greater<>(), eps),
              AlmostEquals(1, 2'532'035));
  EXPECT_THAT(evaluations, Eq(25));

  // 3/4 of the precision.
  eps = Sqrt(Sqrt(ϵ³));
  evaluations = 0;
  EXPECT_THAT(Brent(f, 0.0, π / 2, std::greater<>(), eps),
              AlmostEquals(1, 39'407'194));
  EXPECT_THAT(evaluations, Eq(24));
  evaluations = 0;
  EXPECT_THAT(Brent(f, π / 7, 9 * π / 14, std::greater<>(), eps),
              AlmostEquals(1, 2'528'998));
  EXPECT_THAT(evaluations, Eq(21));

  // 2/3 of the precision.
  eps = Cbrt(ϵ²);
  evaluations = 0;
  EXPECT_THAT(Brent(f, 0.0, π / 2, std::greater<>(), eps),
              AlmostEquals(1, 39'407'194));
  EXPECT_THAT(evaluations, Eq(18));
  evaluations = 0;
  EXPECT_THAT(Brent(f, π / 7, 9 * π / 14, std::greater<>(), eps),
              AlmostEquals(1, 2'461'373));
  EXPECT_THAT(evaluations, Eq(15));

  // 5/9 of the precision.
  eps = Cbrt(Cbrt(ϵ⁵));
  evaluations = 0;
  EXPECT_THAT(Brent(f, 0.0, π / 2, std::greater<>(), eps),
              AlmostEquals(1, 13'620'875));
  EXPECT_THAT(evaluations, Eq(17));
  evaluations = 0;
  EXPECT_THAT(Brent(f, π / 7, 9 * π / 14, std::greater<>(), eps),
              AlmostEquals(1, 18'471'305));
  EXPECT_THAT(evaluations, Eq(11));

  // 1/2 of the precision.
  eps = Sqrt(ϵ);
  evaluations = 0;
  EXPECT_THAT(Brent(f, 0.0, π / 2, std::greater<>(), eps),
              AlmostEquals(1, 44'628'162));
  EXPECT_THAT(evaluations, Eq(9));
  evaluations = 0;
  EXPECT_THAT(Brent(f, π / 7, 9 * π / 14, std::greater<>(), eps),
              AlmostEquals(1, 44'964'716));
  EXPECT_THAT(evaluations, Eq(9));

  // 1/3 of the precision.  We start actually losing precision with respect to
  // the maximum of the theoretical function.
  eps = Cbrt(ϵ);
  evaluations = 0;
  EXPECT_THAT(Brent(f, 0.0, π / 2, std::greater<>(), eps),
              AlmostEquals(1, 8'905'048'021));
  EXPECT_THAT(evaluations, Eq(8));
  evaluations = 0;
  EXPECT_THAT(Brent(f, π / 7, 9 * π / 14, std::greater<>(), eps),
              AlmostEquals(1, 24'970'775));
  EXPECT_THAT(evaluations, Eq(8));
}

TEST_F(RootFindersTest, GoldenSectionSearch) {
  // Arbitrary comparator; we use the lexicographic ordering of the binary
  // representation, with
  // +0 < ... < 1 < ... < +∞ < NaNs < -∞ < ... < -1 < ... < -0.
  EXPECT_THAT(GoldenSectionSearch([](double const x) { return x - 1; },
                                  -π,
                                  π,
                                  [](double const left, double const right) {
                                    return absl::bit_cast<std::uint64_t>(left) <
                                           absl::bit_cast<std::uint64_t>(right);
                                  }),
              AlmostEquals(1, 0));

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
                              t_0 + u * 0.1 * Second,
                              std::less<>()),
          AlmostEquals(t_0 + 3 * π / 2 * Second, 11'863'280, 11'863'284));
    }
  }

  // Maximum.
  EXPECT_THAT(
      GoldenSectionSearch(
          sin,
          t_0 + 1.5 * Second, t_0 + 1.6 * Second, std::greater<>()),
      AlmostEquals(t_0 + π / 2 * Second, 47453133));

  // A big interval will yield a semi-random minimum.
  EXPECT_THAT(GoldenSectionSearch(
                  sin, t_0 - 100 * Second, t_0 + 666 * Second, std::less<>()),
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
  // naïvely, the results are -100000.0031'0988855 and -99999.9968'90111535
  // (separator after the last correct decimal digit), whose error is
  // approximately 3.5 million ULPs.
  auto const s4 = SolveQuadraticEquation(0.0,
                                         1.000'000'000'000'001e25,
                                         2.000'000'000'000'003e20,
                                         1.000'000'000'000'001e15);
  EXPECT_THAT(s4,
              ElementsAre(AlmostEquals(-100000.00316153381194436811, 0),
                          AlmostEquals(-99999.996838466282967631886, 1)));

  // A typed system.
  Instant const t0;
  auto const s5 = SolveQuadraticEquation(
      t0, 1.0 * Metre, 2.0 * Metre / Second, 1.0 * Metre / Second / Second);
  EXPECT_THAT(s5, ElementsAre(t0 - 1.0 * Second));
}

}  // namespace numerics
}  // namespace principia
