#include "numerics/global_optimization.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {

using geometry::Frame;
using geometry::Position;
using quantities::Length;

class GlobalOptimizationTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;
};

TEST_F(GlobalOptimizationTest, Smoke) {
  using Foo = MultiLevelSingleLinkage<Length, Position<World>>;
}

}  // namespace numerics
}  // namespace principia
