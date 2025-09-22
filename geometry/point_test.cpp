#include "geometry/point.hpp"

#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/time_scales.hpp"
#include "base/cpuid.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "numerics/fma.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {

using ::testing::Eq;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_cpuid;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_point;
using namespace principia::geometry::_space;
using namespace principia::numerics::_fma;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

class PointTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  Instant const mjd0 = "MJD0"_TT;
};

using PointDeathTest = PointTest;

TEST_F(PointTest, FMA) {
  if (!CanEmitFMAInstructions || !CPUIDFeatureFlag::FMA.IsSet()) {
    GTEST_SKIP() << "Cannot test FMA on a machine without FMA";
  }
  EXPECT_THAT(FusedMultiplyAdd(3 * Litre, 5 * Second / Litre, mjd0),
              AlmostEquals(mjd0 + 15 * Second, 0));
  EXPECT_THAT(FusedNegatedMultiplyAdd(3, 5 * Second, mjd0),
              AlmostEquals(mjd0 - 15 * Second, 0));
  Displacement<World> const v({-1 * Metre, 2 * Metre, 3 * Metre});
  EXPECT_THAT(FusedMultiplyAdd(2, v, World::origin),
              AlmostEquals(2 * v + World::origin, 0));
  EXPECT_THAT(
      FusedNegatedMultiplyAdd(v / (2 * Second), 3 * Second, World::origin),
      AlmostEquals(World::origin - 1.5 * v, 0));
}

TEST_F(PointTest, Comparisons) {
  EXPECT_TRUE(mjd0 == mjd0);
  EXPECT_FALSE(mjd0 == J2000);
  EXPECT_TRUE(mjd0 != J2000);
  EXPECT_FALSE(mjd0 != mjd0);
}

TEST_F(PointTest, PlusMinus) {
  EXPECT_THAT("MJD0"_TT- "JD0"_TT, Eq(2400000.5 * Day));
  EXPECT_THAT("JD2451545.0"_TT, Eq(J2000));
  EXPECT_THAT("MJD0"_TT - 2400000.5 * Day, Eq("JD0"_TT));
}

TEST_F(PointTest, AssignmentOperators) {
  Instant accumulator = mjd0;
  Instant assignment_result;
  assignment_result = (accumulator += 365 * Day);
  EXPECT_THAT(assignment_result, Eq(accumulator));
  EXPECT_THAT(accumulator, Eq(mjd0 + 365 * Day));
  assignment_result = (accumulator -= 365 * Day);
  EXPECT_THAT(assignment_result, Eq(accumulator));
  EXPECT_THAT(accumulator, Eq(mjd0));
  EXPECT_THAT((accumulator += 365 * Day) -= 365 * Day, Eq(mjd0));
  EXPECT_THAT(accumulator, Eq(mjd0));
}

TEST_F(PointTest, Ordering) {
  // Check that the quantity concept works for double.
  Point<double> zero;
  Point<double> d1 = zero + 1.0;
  Point<double> d2 = zero -3.0;
  EXPECT_TRUE(d2 < d1);
  // Check ordering for instants.
  Instant const t1 = mjd0 + 1 * Day;
  Instant const t2 = mjd0 - 3 * Day;
  EXPECT_TRUE(t2 < t1);
  EXPECT_FALSE(t2 < t2);
  EXPECT_TRUE(t2 <= t1);
  EXPECT_TRUE(t2 <= t2);
  EXPECT_TRUE(t1 > t2);
  EXPECT_FALSE(t1 > t1);
  EXPECT_TRUE(t1 >= t2);
  EXPECT_TRUE(t1 >= t1);
}

// Uncomment to check that non-serializable frames are detected at compile-time.
#if 0
TEST_F(PointTest, SerializationCompilationError) {
  using F = Frame<struct FrameTag>;
  Point<Vector<Length, F>> p;
  serialization::Point message;
  p.WriteToMessage(&message);
}
#endif

TEST_F(PointDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Point message;
    Instant const t0;
    Instant const t1 = t0 + 10 * Second;
    t1.WriteToMessage(&message);
    [[maybe_unused]] Position<World> const d2 =
        Position<World>::ReadFromMessage(message);
  }, "has_multivector");
  EXPECT_DEATH({
    serialization::Point message;
    Position<World> const origin;
    Position<World> const d1 = origin +
        Displacement<World>({-1 * Metre, 2 * Metre, 3 * Metre});
    d1.WriteToMessage(&message);
    [[maybe_unused]] Instant const t2 =
        Instant::ReadFromMessage(message);
  }, "has_scalar");
}

TEST_F(PointTest, SerializationSuccess) {
  serialization::Point message;

  Instant const t0;
  Instant const t1 = t0 + 10 * Second;
  t1.WriteToMessage(&message);
  EXPECT_TRUE(message.has_scalar());
  EXPECT_FALSE(message.has_multivector());
  Instant const t2 = Instant::ReadFromMessage(message);
  EXPECT_EQ(t1, t2);

  Position<World> const origin;
  Position<World> const d1 = origin +
      Displacement<World>({-1 * Metre, 2 * Metre, 3 * Metre});
  d1.WriteToMessage(&message);
  EXPECT_FALSE(message.has_scalar());
  EXPECT_TRUE(message.has_multivector());
  Position<World> const d2 = Position<World>::ReadFromMessage(message);
  EXPECT_EQ(d1, d2);
}

TEST_F(PointDeathTest, BarycentreError) {
  using InstantBarycentreCalculator = BarycentreCalculator<Instant, Volume>;
  EXPECT_DEATH({
    InstantBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
}

TEST_F(PointTest, Barycentres) {
  Instant const t1 = mjd0 + 1 * Day;
  Instant const t2 = mjd0 - 3 * Day;
  Instant const b1 = Barycentre({t1, t2}, {3 * Litre, 1 * Litre});
  Instant const b2 = Barycentre({t2, t1});
  EXPECT_THAT(b1, AlmostEquals(mjd0, 1));
  EXPECT_THAT(b2, Eq(mjd0 - 1 * Day));
}

TEST_F(PointTest, InstantBarycentreCalculator) {
  BarycentreCalculator<Instant, double> calculator;
  Instant const t1 = mjd0 + 2 * Day;
  Instant const t2 = mjd0 - 3 * Day;
  Instant const t3 = mjd0 + 5 * Day;
  Instant const t4 = mjd0 - 7 * Day;
  calculator.Add(t1, 1);
  calculator.Add(t2, 2);
  EXPECT_THAT(calculator.Get(), Eq(mjd0 - 4 * Day / 3));
  calculator.Add(t3, 3);
  calculator.Add(t4, 4);
  EXPECT_THAT(calculator.Get(), Eq(mjd0 - 1.7 * Day));
}

TEST_F(PointTest, DoubleBarycentreCalculator) {
  BarycentreCalculator<Point<double>, double> calculator;
  Point<double> zero;
  Point<double> const d1 = zero + 2;
  Point<double> const d2 = zero - 3;
  Point<double> const d3 = zero + 5;
  Point<double> const d4 = zero - 7;
  calculator.Add(d1, 1);
  calculator.Add(d2, 2);
  EXPECT_THAT(calculator.Get(), Eq(zero - 4.0 / 3.0));
  calculator.Add(d3, 3);
  calculator.Add(d4, 4);
  EXPECT_THAT(calculator.Get(), Eq(zero - 1.7));
}

}  // namespace geometry
}  // namespace principia
