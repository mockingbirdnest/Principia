#include <memory>

#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Geometry/Grassmann.hpp"
#include "Geometry/orthogonal_map.hpp"
#include "Quantities/SI.hpp"

namespace principia {
namespace geometry {

using namespace si;
using namespace testing;

class OrthogonalMapTest : public testing::Test {
 protected:
  struct World;
  typedef OrthogonalMap<World, World> O;

  void SetUp() override {
    vector_.reset(new Vector<quantities::Length, World>(
      R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre)));
    bivector_.reset(new Bivector<quantities::Length, World>(
      R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre)));
    trivector_.reset(new Trivector<quantities::Length, World>(4.0 * Metre));
  }

  std::unique_ptr<Vector<quantities::Length, World>> vector_;
  std::unique_ptr<Bivector<quantities::Length, World>> bivector_;
  std::unique_ptr<Trivector<quantities::Length, World>> trivector_;
};

TEST_F(OrthogonalMapTest, Identity) {
  //EXPECT_THAT(*vector_, MatchesVector(O::Identity()(*vector_), 0.01 * Metre));
}


}  // namespace geometry
}  // namespace principia
