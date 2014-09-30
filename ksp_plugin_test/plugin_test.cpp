
#include "ksp_plugin/plugin.hpp"

#include <memory>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/solar_system.hpp"

using principia::testing_utilities::SolarSystem;

namespace principia {
namespace ksp_plugin {

class PluginTest : public testing::Test {
 protected:
  void SetUp() override {
    solar_system_ = SolarSystem::AtСпутникLaunch();
    bodies_ = solar_system_->massive_bodies();
    initial_time_ = solar_system_->trajectories().front()->last_time();
    sun_index_ = 0;
    sun_gravitational_parameter_ =
      bodies_[sun_index_]->gravitational_parameter();
    planetarium_rotation_ = 0 * SIUnit<Angle>();

    plugin_ = std::make_unique<Plugin>(initial_time_,
                                       sun_index_,
                                       sun_gravitational_parameter_,
                                       planetarium_rotation_);
  }

  std::unique_ptr<SolarSystem> solar_system_;
  SolarSystem::Bodies bodies_;
  Instant initial_time_;
  Index sun_index_;
  GravitationalParameter sun_gravitational_parameter_;
  Angle planetarium_rotation_;

  std::unique_ptr<Plugin> plugin_;
};

TEST_F(PluginTest, DoNothing) {}

}  // namespace ksp_plugin
}  // namespace principia
