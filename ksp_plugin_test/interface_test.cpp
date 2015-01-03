#include "ksp_plugin/interface.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "ksp_plugin/mock_plugin.hpp"

using principia::geometry::Displacement;
using principia::ksp_plugin::AliceSun;
using principia::ksp_plugin::Index;
using principia::ksp_plugin::LineSegment;
using principia::ksp_plugin::MockPlugin;
using principia::ksp_plugin::RenderedTrajectory;
using principia::ksp_plugin::World;
using principia::si::Degree;
using testing::Eq;
using testing::IsNull;
using testing::Ref;
using testing::Return;
using testing::StrictMock;
using testing::_;

bool operator==(XYZ const& left, XYZ const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool operator==(QP const& left, QP const& right) {
  return left.q == right.q && left.p == right.p;
}

namespace {

char const kVesselGUID[] = "NCC-1701-D";

Index const kCelestialIndex = 1;
Index const kParentIndex = 2;

double const kGravitationalParameter = 3;
double const kPlanetariumRotation = 10;
double const kTime = 11;

XYZ kParentPosition = {4, 5, 6};
XYZ kParentVelocity = {7, 8, 9};
QP kParentRelativeDegreesOfFreedom = {kParentPosition, kParentVelocity};

int const kTrajectorySize = 10;

ACTION_TEMPLATE(FillUniquePtr,
                // Note the comma between int and k:
                HAS_1_TEMPLATE_PARAMS(int, k),
                AND_1_VALUE_PARAMS(ptr)) {
  std::tr1::get<k>(args)->reset(ptr);
}

}  // namespace

class InterfaceTest : public testing::Test {
 protected:
  InterfaceTest()
      : plugin_(new StrictMock<MockPlugin>) {}

  std::unique_ptr<StrictMock<MockPlugin>> plugin_;
};

using InterfaceDeathTest = InterfaceTest;

// And there is only one thing we say to Death.
TEST_F(InterfaceDeathTest, Errors) {
  Plugin* plugin = nullptr;
  EXPECT_DEATH({
    principia__DeletePlugin(nullptr);
  }, "pointer.*non NULL");
  EXPECT_DEATH({
    principia__InsertCelestial(plugin,
                               kCelestialIndex,
                               kGravitationalParameter,
                               kParentIndex,
                               kParentRelativeDegreesOfFreedom);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__UpdateCelestialHierarchy(plugin, kCelestialIndex, kParentIndex);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__UpdateCelestialHierarchy(plugin, kCelestialIndex, kParentIndex);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__InsertOrKeepVessel(plugin, kVesselGUID, kParentIndex);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__SetVesselStateOffset(plugin,
                                    kVesselGUID,
                                    kParentRelativeDegreesOfFreedom);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__VesselFromParent(plugin, kVesselGUID);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__CelestialFromParent(plugin, kCelestialIndex);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__NewBodyCentredNonRotatingTransforms(plugin, kCelestialIndex);
  }, "plugin.*non NULL");
}

TEST_F(InterfaceTest, Log) {
  principia__LogInfo("An info");
  principia__LogWarning("A warning");
  principia__LogError("An error");
}

TEST_F(InterfaceTest, NewPlugin) {
  std::unique_ptr<Plugin> plugin(principia__NewPlugin(
                                     kTime,
                                     kParentIndex /*sun_index*/,
                                     kGravitationalParameter,
                                     kPlanetariumRotation));
  EXPECT_THAT(plugin, Not(IsNull()));
}

TEST_F(InterfaceTest, DeletePlugin) {
  Plugin const* plugin = plugin_.release();
  principia__DeletePlugin(&plugin);
  EXPECT_THAT(plugin, IsNull());
}

TEST_F(InterfaceTest, InsertCelestial) {
  EXPECT_CALL(*plugin_,
              InsertCelestial(
                  kCelestialIndex,
                  kGravitationalParameter * SIUnit<GravitationalParameter>(),
                  kParentIndex,
                  RelativeDegreesOfFreedom<AliceSun>(
                      Displacement<AliceSun>(
                          {kParentPosition.x * SIUnit<Length>(),
                           kParentPosition.y * SIUnit<Length>(),
                           kParentPosition.z * SIUnit<Length>()}),
                      Velocity<AliceSun>(
                          {kParentVelocity.x * SIUnit<Speed>(),
                           kParentVelocity.y * SIUnit<Speed>(),
                           kParentVelocity.z * SIUnit<Speed>()}))));
  principia__InsertCelestial(plugin_.get(),
                             kCelestialIndex,
                             kGravitationalParameter,
                             kParentIndex,
                             kParentRelativeDegreesOfFreedom);
}

TEST_F(InterfaceTest, UpdateCelestialHierarchy) {
  EXPECT_CALL(*plugin_,
              UpdateCelestialHierarchy(kCelestialIndex, kParentIndex));
  principia__UpdateCelestialHierarchy(plugin_.get(),
                                      kCelestialIndex,
                                      kParentIndex);
}

TEST_F(InterfaceTest, EndInitialization) {
  EXPECT_CALL(*plugin_,
              EndInitialization());
  principia__EndInitialization(plugin_.get());
}

TEST_F(InterfaceTest, InsertOrKeepVessel) {
  EXPECT_CALL(*plugin_,
              InsertOrKeepVessel(kVesselGUID, kParentIndex));
  principia__InsertOrKeepVessel(plugin_.get(), kVesselGUID, kParentIndex);
}

TEST_F(InterfaceTest, SetVesselStateOffset) {
  EXPECT_CALL(*plugin_,
              SetVesselStateOffset(
                  kVesselGUID,
                  RelativeDegreesOfFreedom<AliceSun>(
                      Displacement<AliceSun>(
                          {kParentPosition.x * SIUnit<Length>(),
                           kParentPosition.y * SIUnit<Length>(),
                           kParentPosition.z * SIUnit<Length>()}),
                      Velocity<AliceSun>(
                          {kParentVelocity.x * SIUnit<Speed>(),
                           kParentVelocity.y * SIUnit<Speed>(),
                           kParentVelocity.z * SIUnit<Speed>()}))));
  principia__SetVesselStateOffset(plugin_.get(),
                                  kVesselGUID,
                                  kParentRelativeDegreesOfFreedom);
}

TEST_F(InterfaceTest, AdvanceTime) {
  EXPECT_CALL(*plugin_,
              AdvanceTime(Instant(kTime * SIUnit<Time>()),
                          kPlanetariumRotation * Degree));
  principia__AdvanceTime(plugin_.get(), kTime, kPlanetariumRotation);
}

TEST_F(InterfaceTest, VesselFromParent) {
  EXPECT_CALL(*plugin_,
              VesselFromParent(kVesselGUID))
      .WillOnce(Return(RelativeDegreesOfFreedom<AliceSun>(
                           Displacement<AliceSun>(
                               {kParentPosition.x * SIUnit<Length>(),
                                kParentPosition.y * SIUnit<Length>(),
                                kParentPosition.z * SIUnit<Length>()}),
                           Velocity<AliceSun>(
                               {kParentVelocity.x * SIUnit<Speed>(),
                                kParentVelocity.y * SIUnit<Speed>(),
                                kParentVelocity.z * SIUnit<Speed>()}))));
  QP const result = principia__VesselFromParent(plugin_.get(),
                                                kVesselGUID);
  EXPECT_THAT(result, Eq(kParentRelativeDegreesOfFreedom));
}

TEST_F(InterfaceTest, CelestialFromParent) {
  EXPECT_CALL(*plugin_,
              CelestialFromParent(kCelestialIndex))
      .WillOnce(Return(RelativeDegreesOfFreedom<AliceSun>(
                           Displacement<AliceSun>(
                               {kParentPosition.x * SIUnit<Length>(),
                                kParentPosition.y * SIUnit<Length>(),
                                kParentPosition.z * SIUnit<Length>()}),
                           Velocity<AliceSun>(
                               {kParentVelocity.x * SIUnit<Speed>(),
                                kParentVelocity.y * SIUnit<Speed>(),
                                kParentVelocity.z * SIUnit<Speed>()}))));
  QP const result = principia__CelestialFromParent(plugin_.get(),
                                                    kCelestialIndex);
  EXPECT_THAT(result, Eq(kParentRelativeDegreesOfFreedom));
}

TEST_F(InterfaceTest, NewBodyCentredNonRotatingFrame) {
  auto dummy_transforms = Transforms<Barycentric, Rendering, Barycentric>::
                              DummyForTesting().release();
  EXPECT_CALL(*plugin_,
              FillBodyCentredNonRotatingTransforms(kCelestialIndex, _))
      .WillOnce(FillUniquePtr<1>(dummy_transforms));
  std::unique_ptr<Transforms<Barycentric, Rendering, Barycentric>> transforms(
      principia__NewBodyCentredNonRotatingTransforms(plugin_.get(),
                                                     kCelestialIndex));
  EXPECT_EQ(dummy_transforms, transforms.get());
}

TEST_F(InterfaceTest, NewBarycentricRotatingFrame) {
  auto dummy_transforms = Transforms<Barycentric, Rendering, Barycentric>::
                              DummyForTesting().release();
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingTransforms(kCelestialIndex,
                                                kParentIndex,
                                                _))
      .WillOnce(FillUniquePtr<2>(dummy_transforms));
  std::unique_ptr<Transforms<Barycentric, Rendering, Barycentric>> transforms(
      principia__NewBarycentricRotatingTransforms(plugin_.get(),
                                                  kCelestialIndex,
                                                  kParentIndex));
  EXPECT_EQ(dummy_transforms, transforms.get());
}

TEST_F(InterfaceTest, DeleteRenderingFrame) {
  auto dummy_transforms = Transforms<Barycentric, Rendering, Barycentric>::
                              DummyForTesting().release();
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingTransforms(kCelestialIndex,
                                                kParentIndex,
                                                _))
      .WillOnce(FillUniquePtr<2>(dummy_transforms));
  Transforms<Barycentric, Rendering, Barycentric>* transforms(
      principia__NewBarycentricRotatingTransforms(plugin_.get(),
                                                  kCelestialIndex,
                                                  kParentIndex));
  EXPECT_EQ(dummy_transforms, transforms);
  principia__DeleteTransforms(&transforms);
  EXPECT_THAT(transforms, IsNull());
}

TEST_F(InterfaceTest, LineAndIterator) {
  auto dummy_transforms = Transforms<Barycentric, Rendering, Barycentric>::
                              DummyForTesting().release();
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingTransforms(kCelestialIndex,
                                                kParentIndex,
                                                _))
      .WillOnce(FillUniquePtr<2>(dummy_transforms));
  Transforms<Barycentric, Rendering, Barycentric>* transforms =
      principia__NewBarycentricRotatingTransforms(plugin_.get(),
                                                  kCelestialIndex,
                                                  kParentIndex);

  // Construct a test rendered trajectory.
  RenderedTrajectory<World> rendered_trajectory;
  Position<World> position =
      World::origin + Displacement<World>(
                          {1 * SIUnit<Length>(),
                           2 * SIUnit<Length>(),
                           3 * SIUnit<Length>()});
  for (int i = 0; i < kTrajectorySize; ++i) {
    Position<World> next_position =
        position + Displacement<World>({10 * SIUnit<Length>(),
                                        20 * SIUnit<Length>(),
                                        30 * SIUnit<Length>()});
    LineSegment<World> line_segment(position, next_position);
    rendered_trajectory.push_back(line_segment);
    position = next_position;
  }

  // Construct a LineAndIterator.
  EXPECT_CALL(*plugin_,
              RenderedVesselTrajectory(
                  kVesselGUID,
                  transforms,
                  World::origin + Displacement<World>(
                                      {kParentPosition.x * SIUnit<Length>(),
                                       kParentPosition.y * SIUnit<Length>(),
                                       kParentPosition.z * SIUnit<Length>()})))
      .WillOnce(Return(rendered_trajectory));
  LineAndIterator* line_and_iterator =
      principia__RenderedVesselTrajectory(plugin_.get(),
                                          kVesselGUID,
                                          transforms,
                                          kParentPosition);
  EXPECT_EQ(kTrajectorySize, line_and_iterator->rendered_trajectory.size());
  EXPECT_EQ(kTrajectorySize, principia__NumberOfSegments(line_and_iterator));

  // Traverse it and check that we get the right data.
  for (int i = 0; i < kTrajectorySize; ++i) {
    EXPECT_FALSE(principia__AtEnd(line_and_iterator));
    XYZSegment const segment = principia__FetchAndIncrement(line_and_iterator);
    EXPECT_EQ(1 + 10 * i, segment.begin.x);
    EXPECT_EQ(2 + 20 * i, segment.begin.y);
    EXPECT_EQ(3 + 30 * i, segment.begin.z);
    EXPECT_EQ(11 + 10 * i, segment.end.x);
    EXPECT_EQ(22 + 20 * i, segment.end.y);
    EXPECT_EQ(33 + 30 * i, segment.end.z);
  }
  EXPECT_TRUE(principia__AtEnd(line_and_iterator));

  // Delete it.
  EXPECT_THAT(line_and_iterator, Not(IsNull()));
  principia__DeleteLineAndIterator(&line_and_iterator);
  EXPECT_THAT(line_and_iterator, IsNull());
}
