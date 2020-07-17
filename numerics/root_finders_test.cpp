
#include "numerics/root_finders.hpp"

#include <vector>

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
  }
}

TEST_F(RootFindersTest, GoldenSectionSearch) {
  Instant const t_0;
  for (int l = 16; l <= 47; ++l) {
    for (int u = 48; u <= 62; ++u) {
      // The result is not overly precise because near its maximum a the
      // function is:
      //   f(a) + fʺ(a) (x - a) / 2 + o((x - a)²)
      // The second order term vanishes when x and a match on the first 26
      // leading bits (roughly).
      EXPECT_THAT(
          GoldenSectionSearch(
              [t_0](Instant const& t) {
                return Sin((t - t_0) * Radian / Second);
              },
              t_0 + l * 0.1 * Second,
              t_0 + u * 0.1 * Second),
          AlmostEquals(t_0 + 3 * π / 2 * Second, 11'863'280, 11'863'284));
    }
  }
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
