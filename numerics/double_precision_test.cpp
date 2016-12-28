
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
using quantities::Length;
using quantities::si::Metre;
using ::testing::Eq;

namespace numerics {
namespace internal_double_precision {

constexpr double ε = std::numeric_limits<double>::epsilon();
constexpr double ε² = ε * ε;

using World = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST,
                    /*inertial=*/false>;

class DoublePrecisionTest : public ::testing::Test {};

TEST_F(DoublePrecisionTest, CompensatedSummation) {
  DoublePrecision<Position<World>> q =
      World::origin + Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre});
  Displacement<World> const δ({ε / 4 * Metre, 0 * Metre, 0 * Metre});
  for (int i = 0; i < 4; ++i) {
    q += δ;
  }
  EXPECT_THAT((q.value - World::origin).coordinates().x, Eq((1 + ε) * Metre));
  EXPECT_THAT(q.error.coordinates().x, Eq(0 * Metre));
}

TEST_F(DoublePrecisionTest, IllConditionedCompensatedSummation) {
  Length const x = (1 + ε) * Metre;
  for (bool cancellation : {true, false}) {
    Length const y = cancellation ? ε * x : π * x;
    DoublePrecision<Length> accumulator;
    accumulator += x;
    accumulator -= y;
    accumulator += x;
    accumulator -= y;
    accumulator -= x;
    accumulator += y;
    accumulator -= x;
    accumulator += y;
    if (cancellation) {
      EXPECT_THAT(accumulator.value, Eq(ε² * Metre));
    } else {
      EXPECT_THAT(accumulator.value, Eq(0 * Metre));
    }
    EXPECT_THAT(accumulator.error, Eq(0 * Metre));
  }
}

TEST_F(DoublePrecisionTest, LongAdd) {
  Length const x = (1 + ε) * Metre;
  for (bool cancellation : {true, false}) {
    Length const y = cancellation ? ε * x : π * x;
    DoublePrecision<Length> accumulator;
    accumulator += TwoSum(x, -y);
    accumulator += TwoSum(x, -y);
    accumulator -= TwoSum(x, -y);
    accumulator -= TwoSum(x, -y);
    if (cancellation) {
      EXPECT_THAT(accumulator.value, Eq(ε² * Metre));
    } else {
      EXPECT_THAT(accumulator.value, Eq(0 * Metre));
    }
    EXPECT_THAT(accumulator.error, Eq(0 * Metre));
  }
}

}  // namespace internal_double_precision
}  // namespace numerics
}  // namespace principia

#include "numerics/double_precision_body.hpp"
