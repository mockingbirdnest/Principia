
#include "numerics/root_finders.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using geometry::Instant;
using quantities::Acceleration;
using quantities::Length;
using quantities::Pow;
using quantities::SIUnit;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Ge;
using ::testing::IsEmpty;
using ::testing::Le;

namespace numerics {

class RootFindersTest : public ::testing::Test {};

// Solving Δt * Δt == n.
TEST_F(RootFindersTest, SquareRoots) {
  Instant const t_0;
  Instant const t_max = t_0 + 10 * Second;
  Length const n_max = Pow<2>(t_max - t_0) * SIUnit<Acceleration>();
  for (Length n = 1 * Metre; n < n_max; n += 1 * Metre) {
    int evaluations = 0;
    auto const equation = [t_0, n, &evaluations](Instant const& t) {
      ++evaluations;
      return Pow<2>(t - t_0) * SIUnit<Acceleration>() - n;
    };
    EXPECT_THAT(Bisect(equation, t_0, t_max) - t_0,
                AlmostEquals(Sqrt(n / SIUnit<Acceleration>()), 0, 1));
    if (n == 25 * Metre) {
      EXPECT_EQ(3, evaluations);
    } else {
      EXPECT_THAT(evaluations, AllOf(Ge(49), Le(58)));
    }
  }
}

TEST_F(RootFindersTest, QuadraticEquations) {
  // Golden ratio.
  auto const s1 = SolveQuadraticEquation(1.0, -1.0, -1.0);
  EXPECT_THAT(s1,
              ElementsAre(AlmostEquals((1 - sqrt(5)) / 2, 1),
                          AlmostEquals((1 + sqrt(5)) / 2, 0)));

  // No solutions.
  auto const s2 = SolveQuadraticEquation(1.0, 0.0, 1.0);
  EXPECT_THAT(s2, IsEmpty());

  // One solution.
  auto const s3 = SolveQuadraticEquation(1.0, 2.0, 1.0);
  EXPECT_THAT(s3, ElementsAre(-1.0));

  // An ill-conditioned system.  I fart in its general direction.
  auto const s4 =
      SolveQuadraticEquation(1.0000001E15, 2.0000003E20, 1.0000001E25);
  EXPECT_THAT(s4,
              ElementsAre(AlmostEquals(-100031.62777541532972762902, 66),
                          AlmostEquals(-99968.38222458367027247098, 65)));
}

}  // namespace numerics
}  // namespace principia
