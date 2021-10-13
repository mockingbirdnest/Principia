#include "physics/discrete_trajectory2.hpp"

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace physics {

using geometry::Frame;

class DiscreteTrajectory2Test : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  DiscreteTrajectory2<World> trajectory_;
};

TEST_F(DiscreteTrajectory2Test) {

}

}  // namespace physics
}  // namespace principia
