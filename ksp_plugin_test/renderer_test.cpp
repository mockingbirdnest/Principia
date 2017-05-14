
#include "ksp_plugin/renderer.hpp"

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin_test/mock_celestial.hpp"
#include "ksp_plugin_test/mock_vessel.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/mock_continuous_trajectory.hpp"
#include "physics/mock_dynamic_frame.hpp"
#include "physics/mock_ephemeris.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_renderer {

using base::not_null;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::MockContinuousTrajectory;
using physics::MockDynamicFrame;
using physics::MockEphemeris;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::_;

class RendererTest : public ::testing::Test {
 protected:
  RendererTest()
      : renderer_(
            &celestial_,
            std::make_unique<MockDynamicFrame<Barycentric, Navigation>>()),
        dynamic_frame_(renderer_.GetPlottingFrame()) {}

  MockCelestial const celestial_;
  Renderer renderer_;
  not_null<MockDynamicFrame<Barycentric, Navigation> const*> const
      dynamic_frame_;
};

TEST_F(RendererTest, TargetVessel) {
  MockEphemeris<Barycentric> ephemeris;
  MockContinuousTrajectory<Barycentric> celestial_trajectory;
  EXPECT_CALL(ephemeris, trajectory(_))
      .WillRepeatedly(Return(&celestial_trajectory));

  MockVessel vessel;
  DiscreteTrajectory<Barycentric> vessel_trajectory;
  vessel_trajectory.Append(Instant{},
                           DegreesOfFreedom<Barycentric>(
                               Barycentric::origin, Velocity<Barycentric>{}));
  EXPECT_CALL(vessel, prediction())
      .WillRepeatedly(ReturnRef(vessel_trajectory));

  renderer_.SetTargetVessel(&vessel, &celestial_, &ephemeris);
  EXPECT_TRUE(renderer_.HasTargetVessel());
  EXPECT_THAT(renderer_.GetTargetVessel(), Ref(vessel));
  MockVessel other_vessel;
  renderer_.ClearTargetVesselIf(&other_vessel);
  EXPECT_THAT(renderer_.GetTargetVessel(), Ref(vessel));
  renderer_.ClearTargetVesselIf(&vessel);
  EXPECT_FALSE(renderer_.HasTargetVessel());
  renderer_.SetTargetVessel(&vessel, &celestial_, &ephemeris);
  renderer_.ClearTargetVessel();
  EXPECT_FALSE(renderer_.HasTargetVessel());
}

}  // namespace internal_renderer
}  // namespace ksp_plugin
}  // namespace principia
