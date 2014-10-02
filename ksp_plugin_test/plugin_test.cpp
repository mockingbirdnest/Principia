
#include "ksp_plugin/plugin.hpp"

#include <map>
#include <memory>

#include "geometry/permutation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/solar_system.hpp"

using principia::geometry::Permutation;
using principia::testing_utilities::AlmostEquals;
using principia::testing_utilities::ICRFJ2000Ecliptic;
using principia::testing_utilities::SolarSystem;
using testing::Eq;

namespace principia {
namespace ksp_plugin {

Index const kSun = 0;
Index const kJupiter = 1;
Index const kSaturn = 2;
Index const kNeptune = 3;
Index const kUranus = 4;
Index const kEarth = 5;
Index const kVenus = 6;
Index const kMars = 7;
Index const kMercury = 8;
Index const kGanymede = 9;
Index const kTitan = 10;
Index const kCallisto = 11;
Index const kIo = 12;
Index const kMoon = 13;
Index const kEuropa = 14;
Index const kTriton = 15;
Index const kEris = 16;
Index const kPluto = 17;

// Body tree.
std::map<Index, Index> const kParents {
    {kJupiter, kSun}, {kSaturn, kSun}, {kNeptune, kSun}, {kUranus, kSun},
    {kEarth, kSun}, {kVenus, kSun}, {kMars, kSun}, {kMercury, kSun},
    {kGanymede, kJupiter}, {kTitan, kSaturn}, {kCallisto, kJupiter},
    {kIo, kJupiter}, {kMoon, kEarth}, {kEuropa, kJupiter}, {kTriton, kNeptune},
    {kEris, kSun}, {kPluto, kSun}};

auto const kLookingGlass = Permutation<ICRFJ2000Ecliptic, AliceSun>(
                               Permutation<ICRFJ2000Ecliptic, AliceSun>::XZY);

class PluginTest : public testing::Test {
 protected:
  void SetUp() override {
    solar_system_ = SolarSystem::AtСпутникLaunch();
    bodies_ = solar_system_->massive_bodies();
    initial_time_ = solar_system_->trajectories().front()->last_time();
    sun_gravitational_parameter_ =
      bodies_[kSun]->gravitational_parameter();
    planetarium_rotation_ = 0 * SIUnit<Angle>();

    plugin_ = std::make_unique<Plugin>(initial_time_,
                                       kSun,
                                       sun_gravitational_parameter_,
                                       planetarium_rotation_);
  }

  std::unique_ptr<SolarSystem> solar_system_;
  SolarSystem::Bodies bodies_;
  Instant initial_time_;
  GravitationalParameter sun_gravitational_parameter_;
  Angle planetarium_rotation_;

  std::unique_ptr<Plugin> plugin_;
};

TEST_F(PluginTest, Initialisation) {
  for (std::size_t index = kSun + 1; index < bodies_.size(); ++index) {
    Index const parent_index = kParents.at(index);
    Displacement<AliceSun> const from_parent_position = kLookingGlass(
        solar_system_->trajectories()[index]->last_position() -
        solar_system_->trajectories()[parent_index]->last_position());
    Velocity<AliceSun> const from_parent_velocity = kLookingGlass(
        solar_system_->trajectories()[index]->last_velocity() -
        solar_system_->trajectories()[parent_index]->last_velocity());
    plugin_->InsertCelestial(index,
                             bodies_[index]->gravitational_parameter(),
                             parent_index,
                             from_parent_position,
                             from_parent_velocity);
  }
  plugin_->EndInitialisation();
  for (std::size_t index = kSun + 1; index < bodies_.size(); ++index) {
    Index const parent_index = kParents.at(index);
    EXPECT_THAT(solar_system_->trajectories()[index]->last_position() -
                solar_system_->trajectories()[parent_index]->last_position(),
                Eq(kLookingGlass.Inverse()(
                       plugin_->CelestialDisplacementFromParent(index))));
    EXPECT_THAT(solar_system_->trajectories()[index]->last_velocity() -
                solar_system_->trajectories()[parent_index]->last_velocity(),
                AlmostEquals(kLookingGlass.Inverse()(
                       plugin_->CelestialParentRelativeVelocity(index))));
  }
}

TEST_F(PluginTest, VesselInsertion) {

}

}  // namespace ksp_plugin
}  // namespace principia
