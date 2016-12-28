
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

}  // namespace internal_double_precision
}  // namespace numerics
}  // namespace principia

#include "numerics/double_precision_body.hpp"
