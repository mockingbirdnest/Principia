
#include "numerics/double_precision.hpp"

#include <limits>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "serialization/geometry.pb.h"

namespace principia {

using geometry::Displacement;
using geometry::Frame;
using geometry::Position;
using geometry::Point;
using quantities::Length;
using quantities::si::Metre;
using ::testing::Eq;

namespace numerics {
namespace internal_double_precision {

constexpr double ε = std::numeric_limits<double>::epsilon();
constexpr double ε² = ε * ε;
constexpr double ε³ = ε² * ε;
constexpr double ε⁴ = ε³ * ε;

using World = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST,
                    /*inertial=*/false>;

class DoublePrecisionTest : public ::testing::Test {};

TEST_F(DoublePrecisionTest, CompensatedSummation) {
  DoublePrecision<Position<World>> accumulator =
      World::origin + Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre});
  Displacement<World> const δ({ε / 4 * Metre, 0 * Metre, 0 * Metre});
  EXPECT_THAT((accumulator + δ).value, Eq(accumulator.value));
  EXPECT_THAT((δ + accumulator).value, Eq(accumulator.value));
  for (int i = 0; i < 4; ++i) {
    accumulator += δ;
  }
  EXPECT_THAT((accumulator.value - World::origin).coordinates().x,
              Eq((1 + ε) * Metre));
  EXPECT_THAT(accumulator.error.coordinates().x, Eq(0 * Metre));
}

TEST_F(DoublePrecisionTest, IllConditionedCompensatedSummation) {
  Length const x = (1 + ε) * Metre;
  Point<Length> const zero;
  for (bool cancellation : {true, false}) {
    Length const y = cancellation ? ε * x : π * x;
    DoublePrecision<Point<Length>> accumulator;
    accumulator += x;
    accumulator -= y;
    accumulator += x;
    accumulator -= y;
    accumulator -= x;
    accumulator += y;
    accumulator -= x;
    accumulator += y;
    if (cancellation) {
      EXPECT_THAT(accumulator.value - zero, Eq(ε² * Metre));
    } else {
      EXPECT_THAT(accumulator.value - zero, Eq(0 * Metre));
    }
    EXPECT_THAT(accumulator.error, Eq(0 * Metre));
  }
}

TEST_F(DoublePrecisionTest, LongAdd) {
  Length const x = (1 + ε) * Metre;
  Point<Length> const zero;
  for (bool cancellation : {true, false}) {
    Length const y = cancellation ? ε * x : π * x;
    DoublePrecision<Point<Length>> accumulator;
    accumulator += TwoSum(+x, -y);
    accumulator -= TwoSum(-x, +y);
    accumulator -= TwoSum(+x, -y);
    accumulator += TwoSum(-x, +y);
    EXPECT_THAT(accumulator.value - zero, Eq(0 * Metre));
    EXPECT_THAT(accumulator.error, Eq(0 * Metre));
  }
}

TEST_F(DoublePrecisionTest, DoubleDoubleDouble) {
  DoublePrecision<DoublePrecision<double>> accumulator;
  accumulator += 1;
  accumulator += ε / 2;
  double const δ = ε² / 4;
  EXPECT_THAT(DebugString(accumulator + δ),
              Eq("+1.00000000000000000e+00|+1.11022302462515654e-16|"
                 "+1.23259516440783095e-32|+0.00000000000000000e+00"));
  for (int i = 0; i < 4; ++i) {
    accumulator += δ;
  }
  accumulator -= ε / 2;
  EXPECT_THAT(DebugString(accumulator),
              Eq("+1.00000000000000000e+00|+4.93038065763132378e-32|"
                 "+0.00000000000000000e+00|+0.00000000000000000e+00"));
  DoublePrecision<DoublePrecision<double>> const long_long_long_float =
      TwoSum(TwoSum(1, ε / 2), TwoSum(ε² / 4, ε³ / 8));
  EXPECT_THAT(DebugString(long_long_long_float),
              Eq("+1.00000000000000000e+00|+1.11022302462515654e-16|"
                 "+1.23259516440783095e-32|+1.36845553156720417e-48"));
  EXPECT_THAT(DebugString(long_long_long_float + TwoSum(0, ε⁴ / 16)),
              Eq(DebugString(long_long_long_float)));
  EXPECT_THAT(DebugString(long_long_long_float + long_long_long_float),
              Eq(DebugString(TwoSum(Scale(2, long_long_long_float.value),
                                    Scale(2, long_long_long_float.error)))));
}

}  // namespace internal_double_precision
}  // namespace numerics
}  // namespace principia

#include "numerics/double_precision_body.hpp"
