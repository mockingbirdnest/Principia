
#include "ksp_plugin/interface.hpp"

#include "base/not_null.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin_test/mock_plugin.hpp"
#include "ksp_plugin_test/mock_renderer.hpp"
#include "ksp_plugin_test/mock_vessel.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/mock_dynamic_frame.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/actions.hpp"

namespace principia {
namespace interface {

using base::make_not_null_unique;
using base::not_null;
using geometry::Rotation;
using ksp_plugin::Barycentric;
using ksp_plugin::Index;
using ksp_plugin::MockPlugin;
using ksp_plugin::MockRenderer;
using ksp_plugin::MockVessel;
using ksp_plugin::Navigation;
using physics::DiscreteTrajectory;
using physics::MockDynamicFrame;
using quantities::si::Metre;
using testing_utilities::FillUniquePtr;
using ::testing::IsNull;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::StrictMock;
using ::testing::_;

namespace {

char const vessel_guid[] = "123-456";

Index const celestial_index = 1;
Index const parent_index = 2;
Index const unused = 666;

int const trajectory_size = 10;

XYZ parent_position = {4, 5, 6};

}  // namespace

class InterfaceRendererTest : public ::testing::Test {
 protected:
  InterfaceRendererTest()
    : plugin_(make_not_null_unique<StrictMock<MockPlugin>>()),
      const_plugin_(plugin_.get()) {}

  not_null<std::unique_ptr<StrictMock<MockPlugin>>> const plugin_;
  StrictMock<MockPlugin> const* const const_plugin_;
  Instant const t0_;
};

TEST_F(InterfaceRendererTest, SetPlottingFrame) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
     mock_navigation_frame =
         new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingNavigationFrame(celestial_index,
                                                     parent_index,
                                                     _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};
  NavigationFrame* navigation_frame =
      principia__NewNavigationFrame(plugin_.get(), parameters);
  EXPECT_EQ(mock_navigation_frame, navigation_frame);
  MockRenderer renderer;
  EXPECT_CALL(*plugin_, renderer()).WillRepeatedly(ReturnRef(renderer));
  EXPECT_CALL(*const_plugin_, renderer()).WillRepeatedly(ReturnRef(renderer));
  EXPECT_CALL(renderer, SetPlottingFrameConstRef(Ref(*navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), &navigation_frame);
  EXPECT_THAT(navigation_frame, IsNull());
  EXPECT_CALL(renderer, GetPlottingFrame())
      .WillOnce(Return(mock_navigation_frame));
  EXPECT_EQ(mock_navigation_frame, principia__GetPlottingFrame(plugin_.get()));
}

TEST_F(InterfaceRendererTest, RenderedPrediction) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
     mock_navigation_frame =
         new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingNavigationFrame(celestial_index,
                                                    parent_index,
                                                    _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};
  NavigationFrame* navigation_frame =
      principia__NewNavigationFrame(plugin_.get(), parameters);
  EXPECT_EQ(mock_navigation_frame, navigation_frame);

  MockRenderer renderer;
  auto const identity = Rotation<Barycentric, AliceSun>::Identity();
  EXPECT_CALL(*plugin_, renderer()).WillRepeatedly(ReturnRef(renderer));
  EXPECT_CALL(*plugin_, PlanetariumRotation())
      .WillRepeatedly(ReturnRef(identity));
  EXPECT_CALL(*plugin_, CurrentTime()).WillOnce(Return(t0_));
  EXPECT_CALL(renderer, SetPlottingFrameConstRef(Ref(*navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), &navigation_frame);
  EXPECT_THAT(navigation_frame, IsNull());

  // Construct a test rendered trajectory.
  auto rendered_trajectory = make_not_null_unique<DiscreteTrajectory<World>>();
  Position<World> position =
      World::origin + Displacement<World>({1 * Metre,
                                           2 * Metre,
                                           3 * Metre});
  rendered_trajectory->Append(
      t0_, DegreesOfFreedom<World>(position, Velocity<World>()));
  for (int i = 1; i < trajectory_size; ++i) {
    position += Displacement<World>({10 * Metre,
                                     20 * Metre,
                                     30 * Metre});
  rendered_trajectory->Append(
      t0_ + i * Second, DegreesOfFreedom<World>(position, Velocity<World>()));
  }

  StrictMock<MockVessel> vessel;
  DiscreteTrajectory<Barycentric> prediction;
  EXPECT_CALL(*plugin_, GetVessel(vessel_guid))
      .WillRepeatedly(Return(&vessel));
  EXPECT_CALL(vessel, prediction())
      .WillRepeatedly(ReturnRef(prediction));
  EXPECT_CALL(renderer,
              FillRenderedBarycentricTrajectoryInWorld(
                  _, _, _,
                  World::origin + Displacement<World>(
                                      {parent_position.x * Metre,
                                       parent_position.y * Metre,
                                       parent_position.z * Metre}),
                  _, _))
      .WillOnce(FillUniquePtr<5>(rendered_trajectory.release()));
  Iterator* iterator =
      principia__RenderedPrediction(plugin_.get(),
                                    vessel_guid,
                                    parent_position);
  EXPECT_EQ(trajectory_size, principia__IteratorSize(iterator));

  // Traverse it and check that we get the right data.
  for (int i = 0; i < trajectory_size; ++i) {
    EXPECT_FALSE(principia__IteratorAtEnd(iterator));
    XYZ const xyz = principia__IteratorGetDiscreteTrajectoryXYZ(iterator);
    EXPECT_EQ(1 + 10 * i, xyz.x);
    EXPECT_EQ(2 + 20 * i, xyz.y);
    EXPECT_EQ(3 + 30 * i, xyz.z);
    principia__IteratorIncrement(iterator);
  }
  EXPECT_TRUE(principia__IteratorAtEnd(iterator));

  // Delete it.
  EXPECT_THAT(iterator, Not(IsNull()));
  principia__IteratorDelete(&iterator);
  EXPECT_THAT(iterator, IsNull());
}

TEST_F(InterfaceRendererTest, Iterator) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
     mock_navigation_frame =
         new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingNavigationFrame(celestial_index,
                                                     parent_index,
                                                     _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};
  NavigationFrame* navigation_frame =
      principia__NewNavigationFrame(plugin_.get(), parameters);
  EXPECT_EQ(mock_navigation_frame, navigation_frame);

  MockRenderer renderer;
  auto const identity = Rotation<Barycentric, AliceSun>::Identity();
  EXPECT_CALL(*plugin_, renderer()).WillRepeatedly(ReturnRef(renderer));
  EXPECT_CALL(*const_plugin_, renderer()).WillRepeatedly(ReturnRef(renderer));
  EXPECT_CALL(*plugin_, PlanetariumRotation())
      .WillRepeatedly(ReturnRef(identity));
  EXPECT_CALL(*plugin_, CurrentTime()).WillOnce(Return(t0_));
  EXPECT_CALL(renderer, SetPlottingFrameConstRef(Ref(*navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), &navigation_frame);
  EXPECT_THAT(navigation_frame, IsNull());

  // Construct a test rendered trajectory.
  auto rendered_trajectory = make_not_null_unique<DiscreteTrajectory<World>>();
  Position<World> position =
      World::origin + Displacement<World>({1 * Metre,
                                           2 * Metre,
                                           3 * Metre});
  rendered_trajectory->Append(
      t0_, DegreesOfFreedom<World>(position, Velocity<World>()));
  for (int i = 1; i < trajectory_size; ++i) {
    position += Displacement<World>({10 * Metre,
                                     20 * Metre,
                                     30 * Metre});
  rendered_trajectory->Append(
      t0_ + i * Second, DegreesOfFreedom<World>(position, Velocity<World>()));
  }

  // Construct a LineAndIterator.
  StrictMock<MockVessel> vessel;
  DiscreteTrajectory<Barycentric> psychohistory;
  EXPECT_CALL(*plugin_, GetVessel(vessel_guid))
      .WillRepeatedly(Return(&vessel));
  EXPECT_CALL(vessel, psychohistory())
      .WillRepeatedly(ReturnRef(psychohistory));
  EXPECT_CALL(renderer,
              FillRenderedBarycentricTrajectoryInWorld(
                  _, _, _,
                  World::origin + Displacement<World>(
                                      {parent_position.x * Metre,
                                       parent_position.y * Metre,
                                       parent_position.z * Metre}),
                  _, _))
      .WillOnce(FillUniquePtr<5>(rendered_trajectory.release()));
  Iterator* iterator =
      principia__RenderedVesselTrajectory(plugin_.get(),
                                          vessel_guid,
                                          parent_position);
  EXPECT_EQ(trajectory_size, principia__IteratorSize(iterator));

  // Traverse it and check that we get the right data.
  for (int i = 0; i < trajectory_size; ++i) {
    EXPECT_FALSE(principia__IteratorAtEnd(iterator));
    XYZ const xyz = principia__IteratorGetDiscreteTrajectoryXYZ(iterator);
    EXPECT_EQ(1 + 10 * i, xyz.x);
    EXPECT_EQ(2 + 20 * i, xyz.y);
    EXPECT_EQ(3 + 30 * i, xyz.z);
    principia__IteratorIncrement(iterator);
  }
  EXPECT_TRUE(principia__IteratorAtEnd(iterator));

  // Delete it.
  EXPECT_THAT(iterator, Not(IsNull()));
  principia__IteratorDelete(&iterator);
  EXPECT_THAT(iterator, IsNull());
}

TEST_F(InterfaceRendererTest, Frenet) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
     mock_navigation_frame =
         new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingNavigationFrame(celestial_index,
                                                     parent_index,
                                                     _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};
  NavigationFrame* navigation_frame =
      principia__NewNavigationFrame(plugin_.get(), parameters);
  EXPECT_EQ(mock_navigation_frame, navigation_frame);

  MockRenderer renderer;
  EXPECT_CALL(*plugin_, renderer()).WillRepeatedly(ReturnRef(renderer));
  EXPECT_CALL(renderer, SetPlottingFrameConstRef(Ref(*navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), &navigation_frame);
  EXPECT_THAT(navigation_frame, IsNull());

  {
    auto const tangent = Vector<double, World>({4, 5, 6});
    EXPECT_CALL(*plugin_, VesselTangent(vessel_guid)).WillOnce(Return(tangent));
    XYZ t = principia__VesselTangent(plugin_.get(), vessel_guid);
    EXPECT_EQ(t.x, tangent.coordinates().x);
    EXPECT_EQ(t.y, tangent.coordinates().y);
    EXPECT_EQ(t.z, tangent.coordinates().z);
  }
  {
    auto const normal = Vector<double, World>({-13, 7, 5});
    EXPECT_CALL(*plugin_, VesselNormal(vessel_guid)).WillOnce(Return(normal));
    XYZ n = principia__VesselNormal(plugin_.get(), vessel_guid);
    EXPECT_EQ(n.x, normal.coordinates().x);
    EXPECT_EQ(n.y, normal.coordinates().y);
    EXPECT_EQ(n.z, normal.coordinates().z);
  }
  {
    auto const binormal = Vector<double, World>({43, 67, 163});
    EXPECT_CALL(*plugin_, VesselBinormal(vessel_guid))
        .WillOnce(Return(binormal));
    XYZ b = principia__VesselBinormal(plugin_.get(), vessel_guid);
    EXPECT_EQ(b.x, binormal.coordinates().x);
    EXPECT_EQ(b.y, binormal.coordinates().y);
    EXPECT_EQ(b.z, binormal.coordinates().z);
  }
  {
    auto const velocity = Velocity<World>(
        {4 * Metre / Second, 5 * Metre / Second, 6 * Metre / Second});
    EXPECT_CALL(*plugin_, VesselVelocity(vessel_guid))
        .WillOnce(Return(velocity));
    XYZ v = principia__VesselVelocity(plugin_.get(), vessel_guid);
    EXPECT_EQ(v.x, velocity.coordinates().x / (Metre / Second));
    EXPECT_EQ(v.y, velocity.coordinates().y / (Metre / Second));
    EXPECT_EQ(v.z, velocity.coordinates().z / (Metre / Second));
  }
}

}  // namespace interface
}  // namespace principia
