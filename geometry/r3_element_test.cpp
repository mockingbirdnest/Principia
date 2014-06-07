#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "quantities/astronomy.hpp"
#include "quantities/constants.hpp"
#include "quantities/bipm.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/uk.hpp"
#include "testing_utilities/algebra.hpp"
#include "testing_utilities/explicit_operators.hpp"

using principia::astronomy::JulianYear;
using principia::astronomy::Parsec;
using principia::bipm::Knot;
using principia::constants::SpeedOfLight;
using principia::quantities::Length;
using principia::quantities::Speed;
using principia::quantities::Time;
using principia::si::Day;
using principia::si::Hour;
using principia::si::Kilo;
using principia::si::Metre;
using principia::si::Minute;
using principia::si::Second;
using principia::testing_utilities::AlmostEquals;
using principia::testing_utilities::Times;
using principia::uk::Furlong;
using principia::uk::Mile;
using principia::uk::Rod;

namespace principia {
namespace geometry {

class R3ElementTest : public testing::Test {
 protected:
  R3Element<Speed> const null_velocity_ = {0 * Knot, 0 * Knot, 0 * Knot};
  R3Element<Speed> const u_ = {3 * Knot, -42 * Parsec / JulianYear, 0 * Knot};
  R3Element<Speed> const v_ = {-π * SpeedOfLight,
                               -e * Kilo(Metre) / Hour,
                               -1 * Knot};
  R3Element<Speed> const w_ = {2 * Mile / Hour,
                               2 * Furlong / Day,
                               2 * Rod / Minute};
  R3Element<Speed> const a_ = {88 * Mile / Hour,
                               300 * Metre / Second,
                               46 * Knot};
};

TEST_F(R3ElementTest, Dumb3Vector) {
  EXPECT_EQ((e * 42) * v_, e * (42 * v_));
  EXPECT_THAT(303.492345479576 * Metre / Second, AlmostEquals(a_.Norm(), 8));
  testing_utilities::TestEquality(42 * v_, 43 * v_);
  testing_utilities::TestVectorSpace<R3Element<Speed>, double>(
      null_velocity_, u_, v_, w_, 0.0, 1.0, e, 42.0, 2);
  testing_utilities::TestAlternatingBilinearMap(
      Cross<Speed, Speed>, u_, v_, w_, a_, 42.0, 2);
  EXPECT_EQ(Cross(R3Element<double>(1, 0, 0),
                  R3Element<double>(0, 1, 0)),
            R3Element<double>(0, 0, 1));
  testing_utilities::TestSymmetricPositiveDefiniteBilinearMap(
      Dot<Speed, Speed>, u_, v_, w_, a_, 42.0, 2);
}

TEST_F(R3ElementTest, MixedProduct) {
  testing_utilities::TestBilinearMap(
      Times<R3Element<Length>, Time, R3Element<Speed>>,
      1 * Second, 1 * JulianYear, u_, v_, 42.0, 2);
  testing_utilities::TestBilinearMap(
      Times<R3Element<Length>, R3Element<Speed>, Time>, w_, a_,
      -1 * Day, 1 * Parsec / SpeedOfLight, -π, 2);
  Time const t = -3 * Second;
  EXPECT_EQ(t * u_, u_ * t);
  EXPECT_THAT((u_ * t) / t, AlmostEquals(u_, 2));
}

}  // namespace geometry
}  // namespace principia
