#include "ksp_plugin/interface.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "ksp_plugin/mock_plugin.hpp"
#include "testing_utilities/death_message.hpp"

using principia::quantities::GravitationalParameter;
using principia::si::Radian;
using principia::si::Second;
using principia::testing_utilities::DeathMessage;
using testing::IsNull;
using testing::StrictMock;

namespace principia {
namespace ksp_plugin {

namespace {

Index const kCelestialIndex = 1;
Index const kParentIndex = 2;

double const kGravitationalParameter = 3;

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

}  // namespace ksp_plugin
}  // namespace principia
