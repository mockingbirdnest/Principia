#include "ksp_plugin/plugin.hpp"

#include <memory>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace ksp_plugin {

class PluginTest : public testing::Test {
 protected:
  void SetUp() override {
    // TODO(phl): Should get these data from SolarSystem.
    initial_time_ = Instant(2 * SIUnit<Time>());
    sun_index_ = 42;
    sun_gravitational_parameter_ = 3 * SIUnit<GravitationalParameter>();
    planetarium_rotation_ = 4 * SIUnit<Angle>();

    plugin_ = std::make_unique<Plugin>(initial_time_,
                                       sun_index_,
                                       sun_gravitational_parameter_,
                                       planetarium_rotation_);
  }

  Instant initial_time_;
  Index sun_index_;
  GravitationalParameter sun_gravitational_parameter_;
  Angle planetarium_rotation_;

  std::unique_ptr<Plugin> plugin_;
};

TEST_F(PluginTest, DoNothing) {}

}  // namespace ksp_plugin
}  // namespace principia
