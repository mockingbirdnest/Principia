
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
using geometry::Displacement;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::MockContinuousTrajectory;
using physics::MockDynamicFrame;
using physics::MockEphemeris;
using physics::RigidMotion;
using physics::RigidTransformation;
using quantities::Time;
using quantities::si::Metre;
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
  void FillTrajectory(
      Instant const& time,
      Time const& step,
      int const number_of_steps,
      std::function<Position<Barycentric>(Instant const& t)> const&
          position_function,
      std::function<Velocity<Barycentric>(Instant const& t)> const&
          velocity_function,
      DiscreteTrajectory<Barycentric>& trajectory) {
    for (int i = 0; i < number_of_steps; ++i) {
      Instant const ti = time + i * step;
      trajectory.Append(ti,
                        DegreesOfFreedom<Barycentric>(position_function(ti),
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
  FillTrajectory(
      /*time=*/t0_,
      /*step=*/1 * Second,
      /*number_of_steps=*/1,
      /*position_function=*/
          [](Instant const& t) { return Barycentric::origin; },
      /*velocity_function=*/
          [](Instant const& t) { return Velocity<Barycentric>(); },
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
  FillTrajectory(
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
      AngularVelocity<Barycentric>(),
      Velocity<Barycentric>());
  for (Instant t = t0_; t < t0_ + 10 * Second; t += 1 * Second) {
    EXPECT_CALL(*dynamic_frame_, ToThisFrameAtTime(t))
        .WillOnce(Return(rigid_motion));
  }

  auto const rendered_trajectory =
      renderer_.RenderBarycentricTrajectoryInPlotting(
          trajectory_to_render.Begin(),
          trajectory_to_render.End());

  EXPECT_EQ(10, rendered_trajectory->Size());
  int index = 0;
  for (auto it = rendered_trajectory->Begin();
       it != rendered_trajectory->End();
       ++it) {
    EXPECT_EQ(Instant{} + index * Second, it.time());
    EXPECT_THAT(
        it.degrees_of_freedom(),
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

  MockVessel vessel;
  DiscreteTrajectory<Barycentric> vessel_trajectory;
  FillTrajectory(
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

  DiscreteTrajectory<Barycentric> trajectory_to_render;
  FillTrajectory(
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
      AngularVelocity<Barycentric>(),
      Velocity<Barycentric>());
  for (Instant t = t0_ + 3 * Second; t < t0_ + 8 * Second; t += 1 * Second) {
    EXPECT_CALL(*dynamic_frame_, ToThisFrameAtTime(t))
        .WillOnce(Return(rigid_motion));
  }

  renderer_.SetTargetVessel(&vessel, &celestial_, &ephemeris);
  auto const rendered_trajectory =
      renderer_.RenderBarycentricTrajectoryInPlotting(
          trajectory_to_render.Begin(),
          trajectory_to_render.End());

  EXPECT_EQ(5, rendered_trajectory->Size());
  int index = 3;
  for (auto it = rendered_trajectory->Begin();
       it != rendered_trajectory->End();
       ++it) {
    EXPECT_EQ(Instant{} + 3 * Second, it.time());
    EXPECT_THAT(
        it.degrees_of_freedom(),
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

}  // namespace internal_renderer
}  // namespace ksp_plugin
}  // namespace principia
