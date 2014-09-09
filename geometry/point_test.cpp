#include "geometry/epoch.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {

using quantities::Time;
using quantities::Volume;
using si::Day;
using si::Litre;
using testing::Eq;
using testing_utilities::AlmostEquals;

class AffineSpaceTest : public testing::Test {
 protected:
};

TEST_F(AffineSpaceTest, Comparisons) {
  EXPECT_TRUE(kUnixEpoch == kUnixEpoch);
  EXPECT_FALSE(kUnixEpoch == kJ2000);
  EXPECT_TRUE(kUnixEpoch != kJ2000);
  EXPECT_FALSE(kUnixEpoch != kUnixEpoch);
}

TEST_F(AffineSpaceTest, PlusMinus) {
  EXPECT_THAT(ModifiedJulianDate(0) - JulianDate(0), Eq(2400000.5 * Day));
  EXPECT_THAT(JulianDate(2451545.0), Eq(kJ2000));
  EXPECT_THAT(ModifiedJulianDate(0) - 2400000.5 * Day, Eq(JulianDate(0)));
}

TEST_F(AffineSpaceTest, AssignmentOperators) {
  Instant accumulator = kUnixEpoch;
  Instant assignment_result;
  assignment_result = (accumulator += 365 * Day);
  EXPECT_THAT(assignment_result, Eq(accumulator));
  EXPECT_THAT(accumulator, Eq(kUnixEpoch + 365 * Day));
  assignment_result = (accumulator -= 365 * Day);
  EXPECT_THAT(assignment_result, Eq(accumulator));
  EXPECT_THAT(accumulator, Eq(kUnixEpoch));
  EXPECT_THAT((accumulator += 365 * Day) -= 365 * Day, Eq(kUnixEpoch));
  EXPECT_THAT(accumulator, Eq(kUnixEpoch));
}

TEST_F(AffineSpaceTest, Ordering) {
  // Check that is_quantity works for double.
  Point<double> d1(1.0);
  Point<double> d2(-3.0);
  EXPECT_TRUE(d2 < d1);
  // Check ordering for instants.
  Instant const t1 = kUnixEpoch + 1 * Day;
  Instant const t2 = kUnixEpoch - 3 * Day;
  EXPECT_TRUE(t2 < t1);
  EXPECT_FALSE(t2 < t2);
  EXPECT_TRUE(t2 <= t1);
  EXPECT_TRUE(t2 <= t2);
  EXPECT_TRUE(t1 > t2);
  EXPECT_FALSE(t1 > t1);
  EXPECT_TRUE(t1 >= t2);
  EXPECT_TRUE(t1 >= t1);
}

TEST_F(AffineSpaceTest, Barycentres) {
  Instant const t1 = kUnixEpoch + 1 * Day;
  Instant const t2 = kUnixEpoch - 3 * Day;
  Instant const b1 = Barycentre<Time, Volume>({t1, t2}, {3 * Litre, 1 * Litre});
  Instant const b2 = Barycentre<Time, double>({t2, t1}, {1, 1});
  EXPECT_THAT(b1, Eq(kUnixEpoch));
  EXPECT_THAT(b2, Eq(kUnixEpoch - 1 * Day));
}

}  // namespace geometry
}  // namespace principia
