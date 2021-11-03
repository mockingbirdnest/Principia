
#include "ksp_plugin/renderer.hpp"

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin_test/mock_celestial.hpp"
#include "ksp_plugin_test/mock_vessel.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_traject0ry.hpp"
#include "physics/mock_continuous_trajectory.hpp"
#include "physics/mock_dynamic_frame.hpp"
#include "physics/mock_ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"

namespace principia {
namespace ksp_plugin {

using base::not_null;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::DefinesFrame;
using geometry::Displacement;
using geometry::Instant;
using geometry::Position;
using geometry::RigidTransformation;
using geometry::Rotation;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::DiscreteTraject0ry;
using physics::MockContinuousTrajectory;
using physics::MockDynamicFrame;
using physics::MockEphemeris;
using physics::RigidMotion;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::AppendTrajectoryTimeline;
using testing_utilities::Componentwise;
using testing_utilities::NewLinearTrajectoryTimeline;
using ::testing::Ref;
using ::testing::Return;
using ::testing::_;

class RendererTest : public ::testing::Test {
 protected:
  RendererTest()
      : renderer_(
            &celestial_,
            std::make_unique<MockDynamicFrame<Barycentric, Navigation>>()),
        dynamic_frame_(renderer_.GetPlottingFrame()) {}

  Instant const t0_;
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
  DiscreteTraject0ry<Barycentric> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(/*v=*/Barycentric::unmoving,
                                  /*Δt=*/1 * Second,
                                  /*t1=*/t0_,
                                  /*t2=*/t0_ + 1 * Second),
      /*to=*/vessel_trajectory);

  EXPECT_CALL(vessel, prediction())
      .WillRepeatedly(Return(vessel_trajectory.segments().begin()));

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

TEST_F(RendererTest, RenderBarycentricTrajectoryInPlottingWithoutTargetVessel) {
  auto const vx = 6 * Metre / Second;
  auto const vy = 5 * Metre / Second;
  auto const vz = 4 * Metre / Second;
  Velocity<Barycentric> const v({vx, vy, vz});
  DiscreteTraject0ry<Barycentric> trajectory_to_render;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(v,
                                  /*Δt=*/1 * Second,
                                  /*t1=*/t0_,
                                  /*t2=*/t0_ + 10 * Second),
      /*to=*/trajectory_to_render);

  RigidMotion<Barycentric, Navigation> rigid_motion(
      RigidTransformation<Barycentric, Navigation>::Identity(),
      Barycentric::nonrotating,
      Barycentric::unmoving);
  for (Instant t = t0_; t < t0_ + 10 * Second; t += 1 * Second) {
    EXPECT_CALL(*dynamic_frame_, ToThisFrameAtTime(t))
        .WillOnce(Return(rigid_motion));
  }

  auto const rendered_trajectory =
      renderer_.RenderBarycentricTrajectoryInPlotting(
          trajectory_to_render.begin(),
          trajectory_to_render.end());

  EXPECT_EQ(10, rendered_trajectory.size());
  int index = 0;
  for (auto const& [time, degrees_of_freedom] : rendered_trajectory) {
    EXPECT_EQ(t0_ + index * Second, time);
    EXPECT_THAT(degrees_of_freedom,
                Componentwise(
                    AlmostEquals(Navigation::origin + Displacement<Navigation>(
                                                          {vx * (time - t0_),
                                                           vy * (time - t0_),
                                                           vz * (time - t0_)}),
                                 0),
                    AlmostEquals(Velocity<Navigation>({vx, vy, vz}), 0)));
    ++index;
  }
}

TEST_F(RendererTest, RenderBarycentricTrajectoryInPlottingWithTargetVessel) {
  MockEphemeris<Barycentric> ephemeris;
  MockContinuousTrajectory<Barycentric> celestial_trajectory;
  EXPECT_CALL(ephemeris, trajectory(_))
      .WillRepeatedly(Return(&celestial_trajectory));

  DiscreteTraject0ry<Barycentric> trajectory_to_render;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          /*v=*/Velocity<Barycentric>(
              {6 * Metre / Second, 5 * Metre / Second, 4 * Metre / Second}),
          /*Δt=*/1 * Second,
          /*t1=*/t0_,
          /*t2=*/t0_ + 10 * Second),
      /*to=*/trajectory_to_render);

  // The prediction is shorter than the |trajectory_to_render|.
  MockVessel vessel;
  DiscreteTraject0ry<Barycentric> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          /*v=*/Velocity<Barycentric>(
              {1 * Metre / Second, 2 * Metre / Second, 3 * Metre / Second}),
          /*Δt=*/1 * Second,
          /*t1=*/t0_ + 3 * Second,
          /*t2=*/t0_ + 8 * Second),
      vessel_trajectory);
  EXPECT_CALL(vessel, prediction())
      .WillRepeatedly(Return(vessel_trajectory.segments().begin()));

  for (Instant t = t0_ + 3 * Second; t < t0_ + 8 * Second; t += 1 * Second) {
    EXPECT_CALL(celestial_trajectory, EvaluateDegreesOfFreedom(t))
        .WillOnce(Return(DegreesOfFreedom<Barycentric>(
            Barycentric::origin + Displacement<Barycentric>(
                                      {300 * Metre, 200 * Metre, 100 * Metre}),
            Barycentric::unmoving)));
  }

  renderer_.SetTargetVessel(&vessel, &celestial_, &ephemeris);
  auto const rendered_trajectory =
      renderer_.RenderBarycentricTrajectoryInPlotting(
          trajectory_to_render.begin(),
          trajectory_to_render.end());

  EXPECT_EQ(5, rendered_trajectory.size());
  int index = 3;
  for (auto const& [time, degrees_of_freedom] : rendered_trajectory) {
    EXPECT_EQ(t0_ + index * Second, time);
    // The degrees of freedom are computed using a real dynamic frame, not a
    // mock.  No point in re-doing the computation here, we just check that the
    // numbers are reasonable.
    EXPECT_LT((degrees_of_freedom.position() - Navigation::origin).Norm(),
              42 * Metre);
    EXPECT_LT(degrees_of_freedom.velocity().Norm(), 6 * Metre / Second);
    ++index;
  }
}

TEST_F(RendererTest, RenderPlottingTrajectoryInWorldWithoutTargetVessel) {
  DiscreteTraject0ry<Navigation> trajectory_to_render;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          /*v=*/Velocity<Navigation>(
              {6 * Metre / Second, 5 * Metre / Second, 4 * Metre / Second}),
          /*Δt=*/1 * Second,
          /*t1=*/t0_,
          /*t2=*/t0_ + 10 * Second),
      trajectory_to_render);

  Instant const rendering_time = t0_ + 5 * Second;
  Position<World> const sun_world_position =
      World::origin +
      Displacement<World>({300 * Metre, 200 * Metre, 100 * Metre});
  Rotation<Barycentric, AliceSun> const planetarium_rotation(
      1 * Radian,
      Bivector<double, Barycentric>({1.0, 1.1, 1.2}),
      DefinesFrame<AliceSun>{});
  RigidMotion<Navigation, Barycentric> rigid_motion(
      RigidTransformation<Navigation, Barycentric>::Identity(),
      Navigation::nonrotating,
      Navigation::unmoving);
  EXPECT_CALL(*dynamic_frame_, FromThisFrameAtTime(rendering_time))
      .WillOnce(Return(rigid_motion));
  EXPECT_CALL(celestial_, current_position(rendering_time))
      .WillOnce(Return(Barycentric::origin));

  auto const rendered_trajectory =
      renderer_.RenderPlottingTrajectoryInWorld(rendering_time,
                                                trajectory_to_render.begin(),
                                                trajectory_to_render.end(),
                                                sun_world_position,
                                                planetarium_rotation);

  EXPECT_EQ(10, rendered_trajectory.size());
  int index = 0;
  for (auto const& [time, degrees_of_freedom] : rendered_trajectory) {
    EXPECT_EQ(t0_ + index * Second, time);
    // The degrees of freedom are computed using real geometrical transforms.
    // No point in re-doing the computation here, we just check that the numbers
    // are reasonable.
    EXPECT_LT((degrees_of_freedom.position() - World::origin).Norm(),
              452 * Metre);
    EXPECT_LT(degrees_of_freedom.velocity().Norm(), 9 * Metre / Second);
    ++index;
  }
}

TEST_F(RendererTest, Serialization) {
  serialization::Renderer message;
  EXPECT_CALL(*dynamic_frame_, WriteToMessage(_));
  renderer_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_plotting_frame());
}

}  // namespace ksp_plugin
}  // namespace principia
