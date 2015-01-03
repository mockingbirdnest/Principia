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

namespace {

char const kVesselGUID[] = "NCC-1701-D";

Index const kCelestialIndex = 1;
Index const kParentIndex = 2;

double const kGravitationalParameter = 3;
double const kPlanetariumRotation = 10;
double const kTime = 11;

XYZ kParentPosition = {4, 5, 6};
XYZ kParentVelocity = {7, 8, 9};

int const kTrajectorySize = 10;

ACTION_P(FillUniquePtr1, p) { return arg1->reset(p); }

ACTION_TEMPLATE(FillUniquePtr,
                // Note the comma between int and k:
                HAS_1_TEMPLATE_PARAMS(int, k),
                AND_1_VALUE_PARAMS(ptr)) {
  std::tr1::get<k>(args)->reset(ptr);
}

class InterfaceTest : public testing::Test {
 protected:
  InterfaceTest()
      : body_centred_non_rotating_frame_placeholder_(nullptr),
        barycentric_rotating_frame_placeholder_(nullptr),
        plugin_(new StrictMock<MockPlugin>) {}

  // Since we cannot create a BodyCentredNonRotatingFrame or a
  // BarycentricRotatingFrame here, we just create a pointer to a bunch of bytes
  // having the right size.
  template<typename T>
  static T* NewPlaceholder() {
    return reinterpret_cast<T*>(new char[sizeof(T)]);
  }

  // These two pointers are not owned.
  BodyCentredNonRotatingFrame* body_centred_non_rotating_frame_placeholder_;
  BarycentricRotatingFrame* barycentric_rotating_frame_placeholder_;

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
                    kParentPosition,
                    kParentVelocity);
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
                                    kParentPosition,
                                    kParentVelocity);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__VesselDisplacementFromParent(plugin, kVesselGUID);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__VesselParentRelativeVelocity(plugin, kVesselGUID);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__CelestialDisplacementFromParent(plugin, kCelestialIndex);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__CelestialParentRelativeVelocity(plugin, kCelestialIndex);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__NewBodyCentredNonRotatingFrame(plugin, kCelestialIndex);
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
                             kParentPosition,
                             kParentVelocity);
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
                                  kParentPosition,
                                  kParentVelocity);
}

TEST_F(InterfaceTest, AdvanceTime) {
  EXPECT_CALL(*plugin_,
              AdvanceTime(Instant(kTime * SIUnit<Time>()),
                          kPlanetariumRotation * Degree));
  principia__AdvanceTime(plugin_.get(), kTime, kPlanetariumRotation);
}

TEST_F(InterfaceTest, VesselDisplacementFromParent) {
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
  XYZ const result = principia__VesselDisplacementFromParent(plugin_.get(),
                                                             kVesselGUID);
  EXPECT_THAT(result, Eq(kParentPosition));
}

TEST_F(InterfaceTest, VesselParentRelativeVelocity) {
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
  XYZ const result = principia__VesselParentRelativeVelocity(plugin_.get(),
                                                             kVesselGUID);
  EXPECT_THAT(result, Eq(kParentVelocity));
}

TEST_F(InterfaceTest, CelestialDisplacementFromParent) {
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
  XYZ const result = principia__CelestialDisplacementFromParent(
                         plugin_.get(),
                         kCelestialIndex);
  EXPECT_THAT(result, Eq(kParentPosition));
}

TEST_F(InterfaceTest, CelestialParentRelativeVelocity) {
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
  XYZ const result = principia__CelestialParentRelativeVelocity(
                         plugin_.get(),
                         kCelestialIndex);
  EXPECT_THAT(result, Eq(kParentVelocity));
}

TEST_F(InterfaceTest, NewBodyCentredNonRotatingFrame) {
  body_centred_non_rotating_frame_placeholder_ =
      NewPlaceholder<BodyCentredNonRotatingFrame>();
  EXPECT_CALL(*plugin_,
              FillBodyCentredNonRotatingFrame(kCelestialIndex, _))
      .WillOnce(FillUniquePtr<1>(body_centred_non_rotating_frame_placeholder_));
  std::unique_ptr<BodyCentredNonRotatingFrame const> frame(
      principia__NewBodyCentredNonRotatingFrame(plugin_.get(),
                                                kCelestialIndex));
  EXPECT_EQ(body_centred_non_rotating_frame_placeholder_, frame.get());
}

TEST_F(InterfaceTest, NewBarycentricRotatingFrame) {
  barycentric_rotating_frame_placeholder_ =
      NewPlaceholder<BarycentricRotatingFrame>();
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingFrame(kCelestialIndex, kParentIndex, _))
      .WillOnce(FillUniquePtr<2>(barycentric_rotating_frame_placeholder_));
  std::unique_ptr<BarycentricRotatingFrame const> frame(
      principia__NewBarycentricRotatingFrame(plugin_.get(),
                                             kCelestialIndex,
                                             kParentIndex));
  EXPECT_EQ(barycentric_rotating_frame_placeholder_, frame.get());
}

TEST_F(InterfaceTest, DeleteRenderingFrame) {
  barycentric_rotating_frame_placeholder_ =
      NewPlaceholder<BarycentricRotatingFrame>();
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingFrame(kCelestialIndex, kParentIndex, _))
      .WillOnce(FillUniquePtr<2>(barycentric_rotating_frame_placeholder_));
  RenderingFrame const* frame =
      principia__NewBarycentricRotatingFrame(plugin_.get(),
                                             kCelestialIndex,
                                             kParentIndex);
  EXPECT_EQ(barycentric_rotating_frame_placeholder_, frame);
  principia__DeleteRenderingFrame(&frame);
  EXPECT_THAT(frame, IsNull());
}

TEST_F(InterfaceTest, LineAndIterator) {
  barycentric_rotating_frame_placeholder_ =
      NewPlaceholder<BarycentricRotatingFrame>();
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingFrame(kCelestialIndex, kParentIndex, _))
      .WillOnce(FillUniquePtr<2>(barycentric_rotating_frame_placeholder_));
  RenderingFrame const* frame =
      principia__NewBarycentricRotatingFrame(plugin_.get(),
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
                  Ref(*frame),
                  World::origin + Displacement<World>(
                                      {kParentPosition.x * SIUnit<Length>(),
                                       kParentPosition.y * SIUnit<Length>(),
                                       kParentPosition.z * SIUnit<Length>()})))
      .WillOnce(Return(rendered_trajectory));
  LineAndIterator* line_and_iterator =
      principia__RenderedVesselTrajectory(plugin_.get(),
                                          kVesselGUID,
                                          frame,
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

}  // namespace
