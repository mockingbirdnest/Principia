#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/algebra.hpp"
#include "testing_utilities/explicit_operators.hpp"

using principia::quantities::Dimensionless;
using principia::quantities::Length;
using principia::quantities::Speed;
using principia::quantities::Time;
using principia::testing_utilities::Times;

namespace principia {
namespace geometry {

class R3ElementTest : public testing::Test {
 protected:
  R3Element<Speed> const null_velocity_ =
      R3Element<Speed>(0 * Knot, 0 * Knot, 0 * Knot);
  R3Element<Speed> const u_ =
      R3Element<Speed>(3 * Knot, -42 * Parsec / JulianYear, 0 * Knot);
  R3Element<Speed> const v_ =
      R3Element<Speed>(-π * SpeedOfLight, -e * Kilo(Metre) / Hour, -1 * Knot);
  R3Element<Speed> const w_ =
      R3Element<Speed>(2 * Mile / Hour, 2 * Furlong / Day, 2 * Rod / Minute);
  R3Element<Speed> const a_ =
      R3Element<Speed>(88 * Mile / Hour, 300 * Metre / Second, 46 * Knot);
};

TEST_F(R3ElementTest, Dumb3Vector) {
  ASSERT_EQ((e * Dimensionless(42)) * v_, e * (Dimensionless(42) * v_));
  ASSERT_EQ(303.492345479576 * Metre / Second, a_.Norm(), 8 * DBL_EPSILON);
  testing_utilities::TestEquality(42 * v_, 43 * v_);
  testing_utilities::TestVectorSpace<R3Element<Speed>, Dimensionless>(
      null_velocity_, u_, v_,
      w_, Dimensionless(0),
      Dimensionless(1), e,
      Dimensionless(42),
      2 * DBL_EPSILON);
  testing_utilities::TestAlternatingBilinearMap(Cross<Speed, Speed>, u_, v_, w_, a_,
                              Dimensionless(42), 2 * DBL_EPSILON);
  ASSERT_EQ(Cross(R3Element<Dimensionless>(1, 0, 0),
                  R3Element<Dimensionless>(0, 1, 0)),
            R3Element<Dimensionless>(0, 0, 1));
  testing_utilities::TestSymmetricPositiveDefiniteBilinearMap(
      Dot<Speed, Speed>, u_, v_, w_, a_, Dimensionless(42), 2 * DBL_EPSILON);
}

TEST_F(R3ElementTest, MixedProduct) {
  testing_utilities::TestBilinearMap(
      Times<R3Element<Length>, Time, R3Element<Speed>>,
      1 * Second, 1 * JulianYear, u_, v_, Dimensionless(42),
      2 * DBL_EPSILON);
  testing_utilities::TestBilinearMap(
      Times<R3Element<Length>, R3Element<Speed>, Time>, w_, a_,
      -1 * Day, 1 * Parsec / SpeedOfLight, Dimensionless(-π),
       2 * DBL_EPSILON);
  Time const t = -3 * Second;
  ASSERT_EQ(t * u_, u_ * t);
  ASSERT_EQ((u_ * t) / t, u_, 2 * DBL_EPSILON);
}

}  // namespace geometry
}  // namespace principia
