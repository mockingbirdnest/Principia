#include "ksp_plugin/interface.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "ksp_plugin/mock_plugin.hpp"
#include "testing_utilities/death_message.hpp"

using principia::si::Degree;
using principia::testing_utilities::DeathMessage;
using testing::Eq;
using testing::IsNull;
using testing::Return;
using testing::StrictMock;

namespace principia {
namespace ksp_plugin {

namespace {

char const kVesselGUID[] = "NCC-1701-D";

Index const kCelestialIndex = 1;
Index const kParentIndex = 2;

double const kGravitationalParameter = 3;
double const kPlanetariumRotation = 10;
double const kTime = 11;

XYZ kParentPosition = {4, 5, 6};
XYZ kParentVelocity = {7, 8, 9};

}  // namespace

class InterfaceTest : public testing::Test {
 protected:
  InterfaceTest()
      : plugin_(new StrictMock<MockPlugin>) {}

  std::unique_ptr<StrictMock<MockPlugin>> plugin_;
};

using InterfaceDeathTest = InterfaceTest;

TEST_F(InterfaceTest, DeletePluginSuccess) {
  Plugin const* plugin = plugin_.release();
  DeletePlugin(&plugin);
  EXPECT_THAT(plugin, IsNull());
}

TEST_F(InterfaceDeathTest, DeletePluginError) {
  EXPECT_DEATH({
    DeletePlugin(nullptr);
  }, DeathMessage("plugin.*non NULL"));
}

TEST_F(InterfaceTest, InsertCelestial) {
  EXPECT_CALL(*plugin_,
              InsertCelestial(
                  kCelestialIndex,
                  kGravitationalParameter * SIUnit<GravitationalParameter>(),
                  kParentIndex,
                  Displacement<AliceSun>(
                      {kParentPosition.x * SIUnit<Length>(),
                       kParentPosition.y * SIUnit<Length>(),
                       kParentPosition.z * SIUnit<Length>()}),
                  Velocity<AliceSun>(
                      {kParentVelocity.x * SIUnit<Speed>(),
                       kParentVelocity.y * SIUnit<Speed>(),
                       kParentVelocity.z * SIUnit<Speed>()})));
  InsertCelestial(plugin_.get(),
                  kCelestialIndex,
                  kGravitationalParameter,
                  kParentIndex,
                  kParentPosition,
                  kParentVelocity);
}

TEST_F(InterfaceTest, UpdateCelestialHierarchy) {
  EXPECT_CALL(*plugin_,
              UpdateCelestialHierarchy(kCelestialIndex, kParentIndex));
  UpdateCelestialHierarchy(plugin_.get(), kCelestialIndex, kParentIndex);
}

TEST_F(InterfaceTest, EndInitialization) {
  EXPECT_CALL(*plugin_,
              EndInitialization());
  EndInitialization(plugin_.get());
}

TEST_F(InterfaceTest, InsertOrKeepVessel) {
  EXPECT_CALL(*plugin_,
              InsertOrKeepVessel(kVesselGUID, kParentIndex));
  InsertOrKeepVessel(plugin_.get(), kVesselGUID, kParentIndex);
}

TEST_F(InterfaceTest, SetVesselStateOffset) {
  EXPECT_CALL(*plugin_,
              SetVesselStateOffset(
                  kVesselGUID,
                  Displacement<AliceSun>(
                      {kParentPosition.x * SIUnit<Length>(),
                       kParentPosition.y * SIUnit<Length>(),
                       kParentPosition.z * SIUnit<Length>()}),
                  Velocity<AliceSun>(
                      {kParentVelocity.x * SIUnit<Speed>(),
                       kParentVelocity.y * SIUnit<Speed>(),
                       kParentVelocity.z * SIUnit<Speed>()})));
  SetVesselStateOffset(plugin_.get(),
                       kVesselGUID,
                       kParentPosition,
                       kParentVelocity);
}

TEST_F(InterfaceTest, AdvanceTime) {
  EXPECT_CALL(*plugin_,
              AdvanceTime(Instant(kTime * SIUnit<Time>()),
                          kPlanetariumRotation * Degree));
  AdvanceTime(plugin_.get(), kTime, kPlanetariumRotation);
}

TEST_F(InterfaceTest, VesselDisplacementFromParent) {
  EXPECT_CALL(*plugin_,
              VesselDisplacementFromParent(kVesselGUID))
      .WillOnce(Return(Displacement<AliceSun>(
                           {kParentPosition.x * SIUnit<Length>(),
                            kParentPosition.y * SIUnit<Length>(),
                            kParentPosition.z * SIUnit<Length>()})));
  XYZ const result = VesselDisplacementFromParent(plugin_.get(), kVesselGUID);
  EXPECT_THAT(result, Eq(kParentPosition));
}

}  // namespace ksp_plugin
}  // namespace principia
