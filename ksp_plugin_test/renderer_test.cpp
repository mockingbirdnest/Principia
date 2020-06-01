
#include "ksp_plugin/renderer.hpp"

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
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
#include "physics/rigid_motion.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_renderer {

using base::not_null;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::DefinesFrame;
using geometry::Displacement;
using geometry::RigidTransformation;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::MockContinuousTrajectory;
using physics::MockDynamicFrame;
using physics::MockEphemeris;
using physics::RigidMotion;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
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

  // TODO(phl): There is a similar function in ContinuousTrajectoryTest and
  // possibly in other places.  It would be good to factor this out.
  template<typename Frame>
  void FillTrajectory(
      Instant const& time,
      Time const& step,
      int const number_of_steps,
      std::function<Position<Frame>(Instant const& t)> const& position_function,
      std::function<Velocity<Frame>(Instant const& t)> const& velocity_function,
      DiscreteTrajectory<Frame>& trajectory) {
    for (int i = 0; i < number_of_steps; ++i) {
      Instant const ti = time + i * step;
      trajectory.Append(ti,
                        DegreesOfFreedom<Frame>(position_function(ti),
                                                velocity_function(ti)));
    }
  }

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
  DiscreteTrajectory<Barycentric> vessel_trajectory;
  FillTrajectory<Barycentric>(
      /*time=*/t0_,
      /*step=*/1 * Second,
      /*number_of_steps=*/1,
      /*position_function=*/
          [](Instant const& t) { return Barycentric::origin; },
      /*velocity_function=*/
          [](Instant const& t) { return Barycentric::unmoving; },
      vessel_trajectory);
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

TEST_F(RendererTest, RenderBarycentricTrajectoryInPlottingWithoutTargetVessel) {
  DiscreteTrajectory<Barycentric> trajectory_to_render;
  FillTrajectory<Barycentric>(
      /*time=*/t0_,
      /*step=*/1 * Second,
      /*number_of_steps=*/10,
      /*position_function=*/
          [this](Instant const& t) {
            return Barycentric::origin +
                   (t - t0_) * Velocity<Barycentric>({6 * Metre / Second,
                                                      5 * Metre / Second,
                                                      4 * Metre / Second});
          },
      /*velocity_function=*/
          [](Instant const& t) {
            return Velocity<Barycentric>(
                {6 * Metre / Second, 5 * Metre / Second, 4 * Metre / Second});
          },
      trajectory_to_render);

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

  EXPECT_EQ(10, rendered_trajectory->Size());
  int index = 0;
  for (auto const& [time, degrees_of_freedom] : *rendered_trajectory) {
    EXPECT_EQ(t0_ + index * Second, time);
    EXPECT_THAT(
        degrees_of_freedom,
        Componentwise(
            AlmostEquals(Navigation::origin +
                             Displacement<Navigation>({6 * index * Metre,
                                                       5 * index * Metre,
                                                       4 * index * Metre}),
                         0),
            AlmostEquals(Velocity<Navigation>({6 * Metre / Second,
                                               5 * Metre / Second,
                                               4 * Metre / Second}),
                         0)));
    ++index;
  }
}

TEST_F(RendererTest, RenderBarycentricTrajectoryInPlottingWithTargetVessel) {
  MockEphemeris<Barycentric> ephemeris;
  MockContinuousTrajectory<Barycentric> celestial_trajectory;
  EXPECT_CALL(ephemeris, trajectory(_))
      .WillRepeatedly(Return(&celestial_trajectory));

  DiscreteTrajectory<Barycentric> trajectory_to_render;
  FillTrajectory<Barycentric>(
      /*time=*/t0_,
      /*step=*/1 * Second,
      /*number_of_steps=*/10,
      /*position_function=*/
          [this](Instant const& t) {
            return Barycentric::origin +
                   (t - t0_) * Velocity<Barycentric>({6 * Metre / Second,
                                                      5 * Metre / Second,
                                                      4 * Metre / Second});
          },
      /*velocity_function=*/
          [](Instant const& t) {
            return Velocity<Barycentric>(
                {6 * Metre / Second, 5 * Metre / Second, 4 * Metre / Second});
          },
      trajectory_to_render);

  // The prediction is shorter than the |trajectory_to_render|.
  MockVessel vessel;
  DiscreteTrajectory<Barycentric> vessel_trajectory;
  FillTrajectory<Barycentric>(
      /*time=*/t0_ + 3 * Second,
      /*step=*/1 * Second,
      /*number_of_steps=*/5,
      /*position_function=*/
      [this](Instant const& t) {
        return Barycentric::origin +
               (t - t0_) * Velocity<Barycentric>({1 * Metre / Second,
                                                  2 * Metre / Second,
                                                  3 * Metre / Second});
      },
      /*velocity_function=*/
      [](Instant const& t) {
        return Velocity<Barycentric>(
            {1 * Metre / Second, 2 * Metre / Second, 3 * Metre / Second});
      },
      vessel_trajectory);
  EXPECT_CALL(vessel, prediction())
      .WillRepeatedly(ReturnRef(vessel_trajectory));

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

  EXPECT_EQ(5, rendered_trajectory->Size());
  int index = 3;
  for (auto const& [time, degrees_of_freedom] : *rendered_trajectory) {
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
  DiscreteTrajectory<Navigation> trajectory_to_render;
  FillTrajectory<Navigation>(
      /*time=*/t0_,
      /*step=*/1 * Second,
      /*number_of_steps=*/10,
      /*position_function=*/
          [this](Instant const& t) {
            return Navigation::origin +
                   (t - t0_) * Velocity<Navigation>({6 * Metre / Second,
                                                     5 * Metre / Second,
                                                     4 * Metre / Second});
          },
      /*velocity_function=*/
          [](Instant const& t) {
            return Velocity<Navigation>(
                {6 * Metre / Second, 5 * Metre / Second, 4 * Metre / Second});
          },
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

  EXPECT_EQ(10, rendered_trajectory->Size());
  int index = 0;
  for (auto const& [time, degrees_of_freedom] : *rendered_trajectory) {
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

}  // namespace internal_renderer
}  // namespace ksp_plugin
}  // namespace principia
