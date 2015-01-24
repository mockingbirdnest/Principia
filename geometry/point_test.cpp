#include "geometry/point.hpp"

#include <vector>

#include "geometry/epoch.hpp"
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
using si::Metre;
using si::Second;
using testing::Eq;
using testing_utilities::AlmostEquals;

class PointTest : public testing::Test {
 protected:
  struct World;
};

using PointDeathTest = PointTest;

TEST_F(PointTest, Comparisons) {
  EXPECT_TRUE(kUnixEpoch == kUnixEpoch);
  EXPECT_FALSE(kUnixEpoch == kJ2000);
  EXPECT_TRUE(kUnixEpoch != kJ2000);
  EXPECT_FALSE(kUnixEpoch != kUnixEpoch);
}

TEST_F(PointTest, PlusMinus) {
  EXPECT_THAT(ModifiedJulianDate(0) - JulianDate(0), Eq(2400000.5 * Day));
  EXPECT_THAT(JulianDate(2451545.0), Eq(kJ2000));
  EXPECT_THAT(ModifiedJulianDate(0) - 2400000.5 * Day, Eq(JulianDate(0)));
}

TEST_F(PointTest, AssignmentOperators) {
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

TEST_F(PointTest, Ordering) {
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

TEST_F(PointDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Point message;
    Instant const t1(10 * Second);
    t1.WriteToMessage(&message);
    Position<World> const d2 = Position<World>::ReadFromMessage(message);
  }, "has_multivector");
  EXPECT_DEATH({
    serialization::Point message;
    Position<World> const d1(
        Displacement<World>({-1 * Metre, 2 * Metre, 3 * Metre}));
    d1.WriteToMessage(&message);
    Instant const t2 = Instant::ReadFromMessage(message);
  }, "has_scalar");
}

TEST_F(PointTest, SerializationSuccess) {
  serialization::Point message;

  Instant const t1(10 * Second);
  t1.WriteToMessage(&message);
  EXPECT_TRUE(message.has_scalar());
  EXPECT_FALSE(message.has_multivector());
  Instant const t2 = Instant::ReadFromMessage(message);
  EXPECT_EQ(t1, t2);

  Position<World> const d1(
      Displacement<World>({-1 * Metre, 2 * Metre, 3 * Metre}));
  d1.WriteToMessage(&message);
  EXPECT_FALSE(message.has_scalar());
  EXPECT_TRUE(message.has_multivector());
  Position<World> const d2 = Position<World>::ReadFromMessage(message);
  EXPECT_EQ(d1, d2);
}

TEST_F(PointDeathTest, BarycentreError) {
  // The <> seem to confuse EXPECT_DEATH, hence the lambda.
  auto barycentre =
      [](std::vector<Instant> const& instants,
         std::vector<Volume> const& weights) -> Instant {
    return Barycentre<Time, Volume>(instants, weights);
  };
  EXPECT_DEATH({
    Instant const t1 = kUnixEpoch + 1 * Day;
    Instant const t2 = kUnixEpoch - 3 * Day;
    barycentre({t1, t2}, {3 * Litre, 4 * Litre, 5 * Litre});
  }, "unequal sizes");
  EXPECT_DEATH({
    barycentre({}, {});
  }, "Empty input");
  EXPECT_DEATH({
    Instant::BarycentreCalculator<Volume> calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
}

TEST_F(PointTest, Barycentres) {
  Instant const t1 = kUnixEpoch + 1 * Day;
  Instant const t2 = kUnixEpoch - 3 * Day;
  Instant const b1 = Barycentre<Time, Volume>({t1, t2}, {3 * Litre, 1 * Litre});
  Instant const b2 = Barycentre<Time, double>({t2, t1}, {1, 1});
  EXPECT_THAT(b1, Eq(kUnixEpoch));
  EXPECT_THAT(b2, Eq(kUnixEpoch - 1 * Day));
}

TEST_F(PointTest, InstantBarycentreCalculator) {
  Instant::BarycentreCalculator<double> calculator;
  Instant const t1 = kUnixEpoch + 2 * Day;
  Instant const t2 = kUnixEpoch - 3 * Day;
  Instant const t3 = kUnixEpoch + 5 * Day;
  Instant const t4 = kUnixEpoch - 7 * Day;
  calculator.Add(t1, 1);
  calculator.Add(t2, 2);
  EXPECT_THAT(calculator.Get(), Eq(kUnixEpoch - 4 * Day / 3));
  calculator.Add(t3, 3);
  calculator.Add(t4, 4);
  EXPECT_THAT(calculator.Get(), Eq(kUnixEpoch - 1.7 * Day));
}

TEST_F(PointTest, DoubleBarycentreCalculator) {
  Point<double>::BarycentreCalculator<double> calculator;
  Point<double> const d1 = Point<double>(2);
  Point<double> const d2 = Point<double>(-3);
  Point<double> const d3 = Point<double>(5);
  Point<double> const d4 = Point<double>(-7);
  calculator.Add(d1, 1);
  calculator.Add(d2, 2);
  EXPECT_THAT(calculator.Get(), Eq(Point<double>(-4.0 / 3.0)));
  calculator.Add(d3, 3);
  calculator.Add(d4, 4);
  EXPECT_THAT(calculator.Get(), Eq(Point<double>(-1.7)));
}

}  // namespace geometry
}  // namespace principia
