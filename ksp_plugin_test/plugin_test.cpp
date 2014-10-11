
#include "ksp_plugin/plugin.hpp"

#include <cmath>
#include <map>
#include <memory>

#include "geometry/permutation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/mock_n_body_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/death_message.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system.hpp"

using principia::geometry::Bivector;
using principia::geometry::Permutation;
using principia::physics::MockNBodySystem;
using principia::quantities::Abs;
using principia::quantities::Sqrt;
using principia::si::Day;
using principia::si::Radian;
using principia::si::AstronomicalUnit;
using principia::testing_utilities::AbsoluteError;
using principia::testing_utilities::RelativeError;
using principia::testing_utilities::AlmostEquals;
using principia::testing_utilities::DeathMessage;
using principia::testing_utilities::ICRFJ2000Ecliptic;
using principia::testing_utilities::SolarSystem;
using testing::AllOf;
using testing::Eq;
using testing::Ge;
using testing::Gt;
using testing::InSequence;
using testing::Lt;
using testing::Ref;
using testing::SizeIs;
using testing::StrictMock;
using testing::_;

namespace principia {
namespace ksp_plugin {

namespace {

// Appends a |DegreesOfFreedom| equal to the last one at the given |time| to
// each |Trajectory| in the |k|th parameter of the expected call.
// |T| must be |NBodySystem<>::Trajectories|, |time| must be an |Instant|.
ACTION_TEMPLATE(AppendTimeToTrajectories,
                HAS_2_TEMPLATE_PARAMS(int, k, typename, T),
                AND_1_VALUE_PARAMS(time)) {
  for (auto* trajectory : T(std::tr1::get<k>(args))) {
    trajectory->Append(time,
                       {trajectory->last_position(),
                        trajectory->last_velocity()});
  }
}

}  // namespace

class TestablePlugin : public Plugin {
 public:
  // Takes ownership of |n_body_system|.
  // TODO(phl): We'd like to pass a |unique_ptr<>| but that seems to confuse
  // gMock.
  TestablePlugin(Instant const& initial_time,
                 Index const sun_index,
                 GravitationalParameter const& sun_gravitational_parameter,
                 Angle const& planetarium_rotation,
                 MockNBodySystem<Barycentre>* n_body_system)
      : Plugin(initial_time,
               sun_index,
               sun_gravitational_parameter,
               planetarium_rotation) {
    n_body_system_.reset(n_body_system);
  }

  Time const& Δt() const {
    return Δt_;
  }

  SPRKIntegrator<Length, Speed> const& prolongation_integrator() const {
    return prolongation_integrator_;
  }

  SPRKIntegrator<Length, Speed> const& history_integrator() const {
    return history_integrator_;
  }
};

class PluginTest : public testing::Test {
 protected:
  PluginTest()
      : looking_glass_(Permutation<ICRFJ2000Ecliptic, AliceSun>::XZY),
        solar_system_(SolarSystem::AtСпутник1Launch(
            SolarSystem::Accuracy::kMajorBodiesOnly)),
        bodies_(solar_system_->massive_bodies()),
        initial_time_(solar_system_->trajectories().front()->last_time()),
        sun_gravitational_parameter_(
            bodies_[SolarSystem::kSun]->gravitational_parameter()),
        planetarium_rotation_(1 * Radian) {
    n_body_system_ = new MockNBodySystem<Barycentre>();
    plugin_ = std::make_unique<StrictMock<TestablePlugin>>(
                  initial_time_,
                  SolarSystem::kSun,
                  sun_gravitational_parameter_,
                  planetarium_rotation_,
                  n_body_system_);
  }

  void InsertAllSolarSystemBodies() {
    for (std::size_t index = SolarSystem::kSun + 1;
         index < bodies_.size();
         ++index) {
      Index const parent_index = SolarSystem::parent(index);
      Displacement<AliceSun> const from_parent_position = looking_glass_(
          solar_system_->trajectories()[index]->last_position() -
          solar_system_->trajectories()[parent_index]->last_position());
      Velocity<AliceSun> const from_parent_velocity = looking_glass_(
          solar_system_->trajectories()[index]->last_velocity() -
          solar_system_->trajectories()[parent_index]->last_velocity());
      plugin_->InsertCelestial(index,
                               bodies_[index]->gravitational_parameter(),
                               parent_index,
                               from_parent_position,
                               from_parent_velocity);
    }
  }

  Permutation<ICRFJ2000Ecliptic, AliceSun> looking_glass_;
  MockNBodySystem<Barycentre>* n_body_system_;  // Not owned.
  std::unique_ptr<SolarSystem> solar_system_;
  SolarSystem::Bodies bodies_;
  Instant initial_time_;
  GravitationalParameter sun_gravitational_parameter_;
  Angle planetarium_rotation_;

  std::unique_ptr<StrictMock<TestablePlugin>> plugin_;
};

TEST_F(PluginTest, Initialization) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    Index const parent_index = SolarSystem::parent(index);
    EXPECT_THAT(solar_system_->trajectories()[index]->last_position() -
                solar_system_->trajectories()[parent_index]->last_position(),
                AlmostEquals(looking_glass_.Inverse()(
                    plugin_->CelestialDisplacementFromParent(index)),
                    250000));
    EXPECT_THAT(solar_system_->trajectories()[index]->last_velocity() -
                solar_system_->trajectories()[parent_index]->last_velocity(),
                AlmostEquals(looking_glass_.Inverse()(
                    plugin_->CelestialParentRelativeVelocity(index)),
                    1000));
  }
}

TEST_F(PluginTest, InsertCelestialError) {
  Displacement<AliceSun> const from_parent_position = looking_glass_(
      solar_system_->trajectories().front()->last_position() -
      solar_system_->trajectories().front()->last_position());
  Velocity<AliceSun> const from_parent_velocity = looking_glass_(
      solar_system_->trajectories().front()->last_velocity() -
      solar_system_->trajectories().front()->last_velocity());
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertCelestial(42,
                             bodies_.front()->gravitational_parameter(),
                             SolarSystem::kSun,
                             from_parent_position,
                             from_parent_velocity);
  }, DeathMessage("before the end of initialization"));
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertCelestial(42,
                             bodies_.front()->gravitational_parameter(),
                             1729,
                             from_parent_position,
                             from_parent_velocity);
  }, DeathMessage("No body at index"));
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertCelestial(SolarSystem::kEarth,
                             bodies_.front()->gravitational_parameter(),
                             SolarSystem::kSun,
                             from_parent_position,
                             from_parent_velocity);
  }, DeathMessage("Body already exists"));
}

TEST_F(PluginTest, VesselInsertionAtInitialization) {
  GUID const guid = "Test Satellite";
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                    SolarSystem::kEarth);
  EXPECT_TRUE(inserted);
  Displacement<AliceSun> const displacement({3111.0 * Kilo(Metre),
                                             4400.0 * Kilo(Metre),
                                             3810.0 * Kilo(Metre)});
  auto const tangent = displacement * Bivector<double, AliceSun>({1, 2, 3});
  Vector<double, AliceSun> unit_tangent = tangent / tangent.Norm();
  EXPECT_THAT(InnerProduct(unit_tangent, displacement / displacement.Norm()),
              Eq(0));
  // This yields a circular orbit.
  Velocity<AliceSun> const velocity =
      Sqrt(bodies_[SolarSystem::kEarth]->gravitational_parameter() /
               displacement.Norm()) * unit_tangent;
  plugin_->SetVesselStateOffset(guid, displacement, velocity);
  EXPECT_THAT(
      AbsoluteError(plugin_->VesselDisplacementFromParent(guid), displacement),
      Lt(DBL_EPSILON * AstronomicalUnit));
  EXPECT_THAT(plugin_->VesselParentRelativeVelocity(guid),
              AlmostEquals(velocity));
}

// Checks that the plugin correctly uses its 10-second-step history even when
// advanced with smaller timesteps.
TEST_F(PluginTest, AdvanceTimeWithCelestials) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  Time const δt = 0.02 * Second;
  Angle const planetarium_rotation = 42 * Radian;
  auto history_step = [this](int const i) {
    return initial_time_ + i * plugin_->Δt(); 
  };
  for (int i = 0; i < 10; ++i) {
    for (Instant t = history_step(i) + δt; t <= history_step(i + 1); t += δt) {
      EXPECT_CALL(*n_body_system_,
                  Integrate(Ref(plugin_->prolongation_integrator()), t,
                            plugin_->Δt(), 0, true, SizeIs(bodies_.size())))
          .RetiresOnSaturation();
      plugin_->AdvanceTime(t, planetarium_rotation);
    }
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->history_integrator()),
                          history_step(i + 1) + δt,
                          plugin_->Δt(), 0, false,
                          SizeIs(bodies_.size())))
        .WillOnce(
             AppendTimeToTrajectories<5, NBodySystem<Barycentre>::Trajectories>(
                 history_step(i + 1)))
        .RetiresOnSaturation();
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->prolongation_integrator()),
                          history_step(i + 1) + δt, plugin_->Δt(), 0, true,
                          SizeIs(bodies_.size())))
        .RetiresOnSaturation();
    plugin_->AdvanceTime(history_step(i + 1) + δt, planetarium_rotation);
  }
}

}  // namespace ksp_plugin
}  // namespace principia
