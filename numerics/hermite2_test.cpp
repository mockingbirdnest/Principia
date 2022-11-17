#include "numerics/hermite2.hpp"

#include <utility>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {

using geometry::Frame;
using geometry::Inertial;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using quantities::Length;
using quantities::si::Metre;
using quantities::si::Second;

namespace numerics {

class Hermite2Test : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag, Inertial>;

  Instant const t0_;
};

TEST_F(Hermite2Test, Precomputed) {
  Hermite2<Length, Instant> h({t0_ + 1 * Second, t0_ + 2 * Second},
                              {33 * Metre, 40 * Metre},
                              -5 * Metre / Second);

  EXPECT_EQ(33 * Metre, h.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(32.5 * Metre, h.Evaluate(t0_ + 1.25 * Second));
  EXPECT_EQ(33.5 * Metre, h.Evaluate(t0_ + 1.5 * Second));
  EXPECT_EQ(36 * Metre, h.Evaluate(t0_ + 1.75 * Second));
  EXPECT_EQ(40 * Metre, h.Evaluate(t0_ + 2 * Second));

  EXPECT_EQ(-5 * Metre / Second, h.EvaluateDerivative(t0_ + 1 * Second));
  EXPECT_EQ(1 * Metre / Second, h.EvaluateDerivative(t0_ + 1.25 * Second));
  EXPECT_EQ(7 * Metre / Second, h.EvaluateDerivative(t0_ + 1.5 * Second));
  EXPECT_EQ(13 * Metre / Second,
            h.EvaluateDerivative(t0_ + 1.75 * Second));
  EXPECT_EQ(19 * Metre / Second, h.EvaluateDerivative(t0_ + 2 * Second));

  EXPECT_EQ(t0_ + 29.0 / 24.0 * Second, h.FindExtremum());
}

TEST_F(Hermite2Test, Typed) {
  // Just here to check that the types work in the presence of affine spaces.
  Hermite2<Position<World>, Instant> h({t0_ + 1 * Second, t0_ + 2 * Second},
                                       {World::origin, World::origin},
                                       World::unmoving);

  EXPECT_EQ(World::origin, h.Evaluate(t0_ + 1.3 * Second));
  EXPECT_EQ(Velocity<World>(), h.EvaluateDerivative(t0_ + 1.7 * Second));
}

}  // namespace numerics
}  // namespace principia
