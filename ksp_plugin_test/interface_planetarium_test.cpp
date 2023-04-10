#include "ksp_plugin/interface.hpp"

#include "geometry/affine_map.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/permutation.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin_test/mock_planetarium.hpp"
#include "ksp_plugin_test/mock_plugin.hpp"
#include "ksp_plugin_test/mock_renderer.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace interface {

using ::testing::ByMove;
using ::testing::IsNull;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::StrictMock;
using ::testing::_;
using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_permutation;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space_transformations;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_planetarium;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin::_renderer;
using namespace principia::quantities::_quantities;

class InterfacePlanetariumTest : public ::testing::Test {
 protected:
  InterfacePlanetariumTest()
      : plugin_(make_not_null_unique<StrictMock<MockPlugin>>()),
        const_plugin_(plugin_.get()) {}

  not_null<std::unique_ptr<StrictMock<MockPlugin>>> const plugin_;
  StrictMock<MockPlugin> const* const const_plugin_;
  Instant const t0_;
};

TEST_F(InterfacePlanetariumTest, ConstructionDestruction) {
  auto const identity = Rotation<Barycentric, AliceSun>::Identity();
  MockRenderer renderer;
  EXPECT_CALL(*const_plugin_, renderer()).WillRepeatedly(ReturnRef(renderer));
  EXPECT_CALL(*plugin_, CurrentTime()).WillRepeatedly(Return(t0_));
  EXPECT_CALL(*plugin_, PlanetariumRotation())
      .WillRepeatedly(ReturnRef(identity));
  EXPECT_CALL(renderer, WorldToPlotting(_, _, _))
      .WillOnce(Return(RigidTransformation<World, Navigation>(
          World::origin,
          Navigation::origin,
          Permutation<World, Navigation>(
              Permutation<World, Navigation>::CoordinatePermutation::YXZ)
              .Forget<OrthogonalMap>()).Forget<Similarity>()));
  EXPECT_CALL(*plugin_, NewPlanetarium(_, _, _))
      .WillOnce(Return(ByMove(std::make_unique<MockPlanetarium>())));

  Planetarium const* planetarium = principia__PlanetariumCreate(plugin_.get(),
                                                                {100, 200, 300},
                                                                {1, 0, 0},
                                                                {0, 1, 0},
                                                                {0, 0, 1},
                                                                {1, 2, 3},
                                                                10,
                                                                90,
                                                                1.0 / 6000,
                                                                {4, 5, 6});
  principia__PlanetariumDelete(&planetarium);
  EXPECT_THAT(planetarium, IsNull());
}

}  // namespace interface
}  // namespace principia
