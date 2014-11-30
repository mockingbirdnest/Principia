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

using R3ElementDeathTest = R3ElementTest;

TEST_F(R3ElementDeathTest, IndexingOperator) {
  EXPECT_DEATH({
    R3Element<Speed> null_velocity = null_velocity_;
    Speed speed = null_velocity[4];
  }, "\\(const int\\)\\:");
  EXPECT_DEATH({
    Speed const& speed = null_velocity_[3];
  }, "\\(const int\\) const\\:");
}

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
  EXPECT_THAT((u_ * t) / t, AlmostEquals(u_, 1));
}

TEST_F(R3ElementDeathTest, OrthogonalizeError) {
  R3Element<Speed> v1 = {1 * Knot, -2 * Knot, 5 * Knot};
  EXPECT_DEATH({
    null_velocity_.Orthogonalize(&v1);
  }, "0.*!= this_norm");
}

TEST_F(R3ElementTest, OrthogonalizeSuccess) {
  R3Element<Length> const v1 = {1 * Metre, -2 * Metre, 5 * Metre};
  R3Element<Length> v2 = {3 * Metre, 4 * Metre, -1 * Metre};
  v1.Orthogonalize(&v2);
  EXPECT_EQ(0 * Metre * Metre, Dot(v1, v2));
  EXPECT_THAT(v2, AlmostEquals(R3Element<Length>({(10.0 / 3.0) * Metre,
                                                  (10.0 / 3.0) * Metre,
                                                  (2.0 / 3.0) * Metre}), 1));
}

TEST_F(R3ElementTest, Normalize) {
  R3Element<Length> const v = {1 * Metre, -2 * Metre, 5 * Metre};
  EXPECT_THAT(Normalize(v),
              AlmostEquals(R3Element<double>({1.0 / sqrt(30.0),
                                              -2.0 / sqrt(30.0),
                                              5.0 / sqrt(30.0)}), 1));
}

}  // namespace geometry
}  // namespace principia
