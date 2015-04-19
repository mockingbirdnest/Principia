
#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "geometry/permutation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/mock_n_body_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system.hpp"

namespace principia {

using geometry::Bivector;
using geometry::Permutation;
using physics::MockNBodySystem;
using quantities::Abs;
using quantities::ArcTan;
using quantities::Cos;
using quantities::Sin;
using quantities::Sqrt;
using si::Day;
using si::Hour;
using si::Minute;
using si::Radian;
using si::AstronomicalUnit;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::ICRFJ2000Ecliptic;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystem;
using ::testing::AllOf;
using ::testing::Contains;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::InSequence;
using ::testing::Le;
using ::testing::Lt;
using ::testing::Ref;
using ::testing::SizeIs;
using ::testing::StrictMock;
using ::testing::_;

namespace ksp_plugin {

namespace {

int const kNotABody = 1729;

// Appends a |DegreesOfFreedom| equal to the last one at the given |time| to
// each |Trajectory| in the |k|th parameter of the expected call.
// This parameter must be an |NBodySystem<Barycentric>::Trajectories|, |time|
// must be an |Instant|.
ACTION_TEMPLATE(AppendTimeToTrajectories,
                HAS_1_TEMPLATE_PARAMS(int, k),
                AND_1_VALUE_PARAMS(time)) {
  for (auto trajectory : static_cast<NBodySystem<Barycentric>::Trajectories>(
                             std::tr1::get<k>(args))) {
    trajectory->Append(time, trajectory->last().degrees_of_freedom());
  }
}

MATCHER_P(HasNonvanishingIntrinsicAccelerationAt, t, "") {
  if (arg->has_intrinsic_acceleration()) {
    if (arg->evaluate_intrinsic_acceleration(t) ==
        Vector<Acceleration, Barycentric>()) {
       *result_listener << "has vanishing intrinsic acceleration at";
       return false;
    } else {
       *result_listener << "has nonvanishing intrinsic acceleration at";
       return true;
    }
  }
  *result_listener << "has no intrinsic acceleration";
  return false;
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
                 not_null<MockNBodySystem<Barycentric>*> const n_body_system)
      : Plugin(initial_time,
               sun_index,
               sun_gravitational_parameter,
               planetarium_rotation) {
    n_body_system_.reset(n_body_system);
  }

  Time const& Δt() const {
    return Δt_;
  }

  SRKNIntegrator const& prolongation_integrator() const {
    return *prolongation_integrator_;
  }

  SRKNIntegrator const& history_integrator() const {
    return *history_integrator_;
  }
};

class PluginTest : public testing::Test {
 protected:
  PluginTest()
      : looking_glass_(Permutation<ICRFJ2000Ecliptic, AliceSun>::XZY),
        n_body_system_(new MockNBodySystem<Barycentric>()),
        solar_system_(SolarSystem::AtСпутник1Launch(
            SolarSystem::Accuracy::kMajorBodiesOnly)),
        bodies_(solar_system_->massive_bodies()),
        initial_time_(solar_system_->trajectories().front()->last().time()),
        sun_gravitational_parameter_(
            bodies_[SolarSystem::kSun]->gravitational_parameter()),
        planetarium_rotation_(1 * Radian),
        plugin_(make_not_null_unique<StrictMock<TestablePlugin>>(
                    initial_time_,
                    SolarSystem::kSun,
                    sun_gravitational_parameter_,
                    planetarium_rotation_,
                    n_body_system_)) {
    satellite_initial_displacement_ =
        Displacement<AliceSun>({3111.0 * Kilo(Metre),
                                4400.0 * Kilo(Metre),
                                3810.0 * Kilo(Metre)});
    auto const tangent =
        satellite_initial_displacement_ * Bivector<double, AliceSun>({1, 2, 3});
    Vector<double, AliceSun> const unit_tangent = Normalize(tangent);
    EXPECT_THAT(
        InnerProduct(unit_tangent,
                     satellite_initial_displacement_ /
                         satellite_initial_displacement_.Norm()),
        Eq(0));
    // This yields a circular orbit.
    satellite_initial_velocity_ =
        Sqrt(bodies_[SolarSystem::kEarth]->gravitational_parameter() /
                 satellite_initial_displacement_.Norm()) * unit_tangent;
  }

  void InsertAllSolarSystemBodies() {
    for (std::size_t index = SolarSystem::kSun + 1;
         index < bodies_.size();
         ++index) {
      Index const parent_index = SolarSystem::parent(index);
      RelativeDegreesOfFreedom<AliceSun> const from_parent = looking_glass_(
          solar_system_->trajectories()[index]->
              last().degrees_of_freedom() -
          solar_system_->trajectories()[parent_index]->
              last().degrees_of_freedom());
      plugin_->InsertCelestial(index,
                               bodies_[index]->gravitational_parameter(),
                               parent_index,
                               from_parent);
    }
  }

  // The time of the |step|th synchronized history step of |plugin_|.
  // |HistoryTime(0)| is |initial_time_|.
  Instant HistoryTime(int const step) {
    return initial_time_ + step * plugin_->Δt();
  }

  // Keeps the vessel with the given |guid| during the next call to
  // |AdvanceTime|.  The vessel must be present.
  void KeepVessel(GUID const& guid) {
    bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                      SolarSystem::kEarth);
    EXPECT_FALSE(inserted) << guid;
  }

  // Inserts a vessel with the given |guid| and makes it a satellite of Earth
  // with relative position |satellite_initial_displacement_| and velocity
  // |satellite_initial_velocity_|.  The vessel must not already be present.
  // Increments the counter |*number_of_new_vessels|.  |number_of_new_vessels|
  // must not be null.
  void InsertVessel(GUID const& guid,
                    std::size_t* const number_of_new_vessels) {
    bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                      SolarSystem::kEarth);
    EXPECT_TRUE(inserted) << guid;
    plugin_->SetVesselStateOffset(guid,
                                  RelativeDegreesOfFreedom<AliceSun>(
                                      satellite_initial_displacement_,
                                      satellite_initial_velocity_));
    ++*CHECK_NOTNULL(number_of_new_vessels);
  }

  Permutation<ICRFJ2000Ecliptic, AliceSun> looking_glass_;
  not_null<MockNBodySystem<Barycentric>*> n_body_system_;  // Not owned.
  not_null<std::unique_ptr<SolarSystem>> solar_system_;
  SolarSystem::Bodies bodies_;
  Instant initial_time_;
  GravitationalParameter sun_gravitational_parameter_;
  Angle planetarium_rotation_;

  not_null<std::unique_ptr<StrictMock<TestablePlugin>>> plugin_;

  // These initial conditions will yield a low circular orbit around Earth.
  Displacement<AliceSun> satellite_initial_displacement_;
  Velocity<AliceSun> satellite_initial_velocity_;
};

using PluginDeathTest = PluginTest;

TEST_F(PluginDeathTest, SerializationError) {
  EXPECT_DEATH({
    auto plugin =
        make_not_null_unique<Plugin>(
            initial_time_,
            SolarSystem::kSun,
            sun_gravitational_parameter_,
            planetarium_rotation_);
    serialization::Plugin message;
    plugin->WriteToMessage(&message);
  }, "!initializing");
}

TEST_F(PluginTest, Serialization) {
  GUID const satellite = "satellite";
  // We need an actual |Plugin| here rather than a |TestablePlugin|, since
  // that's what |ReadFromMessage| returns.
  auto plugin = make_not_null_unique<Plugin>(
                    initial_time_,
                    SolarSystem::kSun,
                    sun_gravitational_parameter_,
                    planetarium_rotation_);
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    Index const parent_index = SolarSystem::parent(index);
    RelativeDegreesOfFreedom<AliceSun> const from_parent= looking_glass_(
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom() -
        solar_system_->trajectories()[parent_index]->
            last().degrees_of_freedom());
    plugin->InsertCelestial(index,
                            bodies_[index]->gravitational_parameter(),
                            parent_index,
                            from_parent);
  }
  plugin->EndInitialization();
  plugin->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  plugin->SetVesselStateOffset(satellite,
                               RelativeDegreesOfFreedom<AliceSun>(
                                   satellite_initial_displacement_,
                                   satellite_initial_velocity_));
  serialization::Plugin message;
  plugin->WriteToMessage(&message);
  plugin = Plugin::ReadFromMessage(message);
  serialization::Plugin second_message;
  plugin->WriteToMessage(&second_message);
  EXPECT_EQ(message.SerializeAsString(), second_message.SerializeAsString());
  EXPECT_EQ(bodies_.size(), message.celestial_size());
  EXPECT_EQ(1, message.vessel_size());
  EXPECT_EQ(SolarSystem::kEarth, message.vessel(0).parent_index());
  EXPECT_FALSE(message.bubble().has_current());
}

TEST_F(PluginTest, Initialization) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    Index const parent_index = SolarSystem::parent(index);
    RelativeDegreesOfFreedom<ICRFJ2000Ecliptic> const from_parent =
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom() -
        solar_system_->trajectories()[parent_index]->
            last().degrees_of_freedom();
    EXPECT_THAT(from_parent,
                Componentwise(
                    AlmostEquals(looking_glass_.Inverse()(
                            plugin_->CelestialFromParent(index).displacement()),
                        1, 216320),
                    AlmostEquals(looking_glass_.Inverse()(
                            plugin_->CelestialFromParent(index).velocity()),
                        0, 936)));
  }
}

TEST_F(PluginDeathTest, InsertCelestialError) {
  RelativeDegreesOfFreedom<AliceSun> const from_parent = looking_glass_(
      solar_system_->trajectories().front()->last().degrees_of_freedom() -
      solar_system_->trajectories().front()->last().degrees_of_freedom());
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertCelestial(42,
                             bodies_.front()->gravitational_parameter(),
                             SolarSystem::kSun,
                             from_parent);
  }, "before the end of initialization");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertCelestial(42,
                             bodies_.front()->gravitational_parameter(),
                             kNotABody,
                             from_parent);
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertCelestial(SolarSystem::kEarth,
                             bodies_.front()->gravitational_parameter(),
                             SolarSystem::kSun,
                             from_parent);
  }, "Body already exists");
}

TEST_F(PluginDeathTest, UpdateCelestialHierarchyError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->UpdateCelestialHierarchy(SolarSystem::kSun, SolarSystem::kPluto);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->UpdateCelestialHierarchy(kNotABody, SolarSystem::kPluto);
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->UpdateCelestialHierarchy(SolarSystem::kSun, kNotABody);
  }, "Map key not found");
}

TEST_F(PluginDeathTest, InsertOrKeepVesselError) {
  GUID const guid = "Syrio Forel";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertOrKeepVessel(guid, SolarSystem::kSun);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, kNotABody);
  }, "Map key not found");
}

TEST_F(PluginDeathTest, SetVesselStateOffsetError) {
  GUID const guid = "Test Satellite";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->SetVesselStateOffset(guid,
                                  RelativeDegreesOfFreedom<AliceSun>(
                                      satellite_initial_displacement_,
                                      satellite_initial_velocity_));
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->SetVesselStateOffset(guid,
                                  RelativeDegreesOfFreedom<AliceSun>(
                                      satellite_initial_displacement_,
                                      satellite_initial_velocity_));
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, SolarSystem::kSun);
    plugin_->SetVesselStateOffset(guid,
                                  RelativeDegreesOfFreedom<AliceSun>(
                                      satellite_initial_displacement_,
                                      satellite_initial_velocity_));
    plugin_->SetVesselStateOffset(guid,
                                  RelativeDegreesOfFreedom<AliceSun>(
                                      satellite_initial_displacement_,
                                      satellite_initial_velocity_));
  }, "already has a trajectory");
}

TEST_F(PluginDeathTest, AdvanceTimeError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->AdvanceTime(Instant(), Angle());
  }, "Check failed: !initializing");
}

TEST_F(PluginDeathTest, VesselFromParentError) {
  GUID const guid = "Test Satellite";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->VesselFromParent(guid);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->VesselFromParent(guid);
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, SolarSystem::kSun);
    plugin_->VesselFromParent(guid);
  }, "not given an initial state");
}

TEST_F(PluginDeathTest, CelestialFromParentError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->CelestialFromParent(SolarSystem::kEarth);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->CelestialFromParent(kNotABody);
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->CelestialFromParent(SolarSystem::kSun);
  }, "is the sun");
}

TEST_F(PluginTest, VesselInsertionAtInitialization) {
  GUID const guid = "Test Satellite";
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                    SolarSystem::kEarth);
  EXPECT_TRUE(inserted);
  plugin_->SetVesselStateOffset(guid,
                                RelativeDegreesOfFreedom<AliceSun>(
                                    satellite_initial_displacement_,
                                    satellite_initial_velocity_));
  EXPECT_THAT(plugin_->VesselFromParent(guid),
              Componentwise(
                  AlmostEquals(satellite_initial_displacement_, 7460),
                  AlmostEquals(satellite_initial_velocity_, 3)));
}

// Checks that the plugin correctly uses its 10-second-step history even when
// advanced with smaller timesteps.
TEST_F(PluginTest, AdvanceTimeWithCelestialsOnly) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  Time const δt = 0.02 * Second;
  Angle const planetarium_rotation = 42 * Radian;
  for (int step = 0; step < 10; ++step) {
    for (Instant t = HistoryTime(step) + 2 * δt;
         t <= HistoryTime(step + 1);
         t += δt) {
      // Called to compute the prolongations.
      EXPECT_CALL(*n_body_system_,
                  Integrate(Ref(plugin_->prolongation_integrator()), t,
                            plugin_->Δt(), 0, true, SizeIs(bodies_.size())))
          .RetiresOnSaturation();
      plugin_->AdvanceTime(t, planetarium_rotation);
    }
    // Called to advance the synchronized histories.
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->history_integrator()),
                          HistoryTime(step + 1) + δt,
                          plugin_->Δt(), 0, false,
                          SizeIs(bodies_.size())))
        .WillOnce(AppendTimeToTrajectories<5>(HistoryTime(step + 1)))
        .RetiresOnSaturation();
    // Called to compute the prolongations.
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->prolongation_integrator()),
                          HistoryTime(step + 1) + δt, plugin_->Δt(), 0, true,
                          SizeIs(bodies_.size())))
        .RetiresOnSaturation();
    plugin_->AdvanceTime(HistoryTime(step + 1) + δt, planetarium_rotation);
  }
}

// Checks that the plugin correctly advances the history of newly inserted
// vessels with the prolongation integrator (using small steps), then switches
// to the history integrator.  Also tests the removal of vessels.
// The sequence of additions and removals is as follows (using the constants
// defined below; recall that |initial_time_ == HistoryTime(0)|).
// * |HistoryTime(0)|               : Insert |enterprise|.
// * |HistoryTime(0) + a_while|     : Insert |enterprise_d|.
// * |HistoryTime(1) + a_while|     : Insert |stargazer| and |bradbury|.
// * |HistoryTime(1) + half_a_step| : Remove |bradbury| (it is present at
//                                    |HistoryTime(1) + half_a_step|, but not at
//                                    |HistoryTime(1) + half_a_step + δt|).
//                                    It is thus removed while unsynchronized.
// * |HistoryTime(2)|               : Insert |constantinople|.
// * |HistoryTime(3)|               : Remove |enterprise|.
TEST_F(PluginTest, AdvanceTimeWithVessels) {
  Time const δt = 0.02 * Second;
  Time const a_while = 10 * δt;
  Time const half_a_step = 250 * δt;
  Time const ε_δt = 0.1 * δt;
  EXPECT_THAT(half_a_step, Eq(plugin_->Δt() / 2));
  GUID const enterprise = "NCC-1701";
  GUID const enterprise_d = "NCC-1701-D";
  GUID const stargazer = "NCC-2893";
  GUID const bradbury = "NX-72307";
  GUID const constantinople = "NCC-43622";
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  Angle const planetarium_rotation = 42 * Radian;
  std::size_t expected_number_of_new_vessels = 0;
  std::size_t expected_number_of_old_vessels = 0;
  InsertVessel(enterprise, &expected_number_of_new_vessels);
  for (int step = 0; step < 10; ++step) {
    for (Instant t = HistoryTime(step) + 2 * δt;
         t <= HistoryTime(step + 1);
         t += δt) {
      // Keep our vessels.  Make sure we're not inserting new ones.
      if (step <= 3) {
        KeepVessel(enterprise);
      }
      if (t > HistoryTime(0) + a_while + ε_δt) {
        KeepVessel(enterprise_d);
      }
      if (t > HistoryTime(1) + a_while + ε_δt) {
        KeepVessel(stargazer);
      }
      if (t > HistoryTime(1) + a_while + ε_δt &&
          t < HistoryTime(1) + half_a_step - ε_δt) {
        KeepVessel(bradbury);
      } else if (AbsoluteError(t - HistoryTime(1), half_a_step) < ε_δt) {
        // We will be removing |bradbury| in this step.
        --expected_number_of_new_vessels;
      }
      if (step > 2) {
        KeepVessel(constantinople);
      }
      // Called to compute the prolongations and advance the unsynchronized
      // histories.
      EXPECT_CALL(*n_body_system_,
                  Integrate(Ref(plugin_->prolongation_integrator()), t,
                            plugin_->Δt(), 0, true,
                            SizeIs(bodies_.size() +
                                       expected_number_of_old_vessels +
                                       expected_number_of_new_vessels)))
          .RetiresOnSaturation();
      plugin_->AdvanceTime(t, planetarium_rotation);
      if (AbsoluteError(t - HistoryTime(0), a_while) < ε_δt) {
        InsertVessel(enterprise_d, &expected_number_of_new_vessels);
      } else if (AbsoluteError(t - HistoryTime(1), a_while) < ε_δt) {
        InsertVessel(stargazer, &expected_number_of_new_vessels);
        InsertVessel(bradbury, &expected_number_of_new_vessels);
      }
    }
    // Keep the vessels for the history-advancing step.
    if (step <= 3) {
      KeepVessel(enterprise);
    }
    KeepVessel(enterprise_d);
    if (step >= 1) {
      KeepVessel(stargazer);
    }
    if (step > 2) {
      KeepVessel(constantinople);
    }
    // Called to advance the synchronized histories.
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->history_integrator()),
                          HistoryTime(step + 1) + δt,
                          plugin_->Δt(), 0, false,
                          SizeIs(bodies_.size() +
                                     expected_number_of_old_vessels)))
        .WillOnce(AppendTimeToTrajectories<5>(HistoryTime(step + 1)))
        .RetiresOnSaturation();
    if (expected_number_of_new_vessels > 0) {
      // Called to synchronize the new histories.
      EXPECT_CALL(*n_body_system_,
                  Integrate(Ref(plugin_->prolongation_integrator()),
                            HistoryTime(step + 1),
                            plugin_->Δt(), 0, true,
                            SizeIs(bodies_.size() +
                                       expected_number_of_new_vessels)))
          .WillOnce(AppendTimeToTrajectories<5>(HistoryTime(step + 1)))
          .RetiresOnSaturation();
    }
    expected_number_of_old_vessels += expected_number_of_new_vessels;
    expected_number_of_new_vessels = 0;
    // Called to compute the prolongations.
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->prolongation_integrator()),
                          HistoryTime(step + 1) + δt, plugin_->Δt(), 0, true,
                          SizeIs(bodies_.size() +
                                     expected_number_of_old_vessels)))
        .RetiresOnSaturation();
    plugin_->AdvanceTime(HistoryTime(step + 1) + δt, planetarium_rotation);
    if (step == 2) {
      InsertVessel(constantinople, &expected_number_of_new_vessels);
    } else if (step == 3) {
      // We will be removing |enterprise|.
      --expected_number_of_old_vessels;
    }
  }
}

// Checks that the plugin correctly advances the history of newly inserted
// vessels with the prolongation integrator (using small steps), then switches
// to the history integrator.  Also tests the removal of vessels.
// The sequence of additions and removals is as follows (using the constants
// defined below; recall that |initial_time_ == HistoryTime(0)|).
// * |HistoryTime(0)|               : Insert |enterprise|, insert into bubble
//                                    (1 part).
// * |HistoryTime(0) + a_while|     : Insert |enterprise_d|.
// * |HistoryTime(1) + a_while|     : Insert |enterprise_d| into bubble
//                                    (2 parts).
// * |HistoryTime(1) + half_a_step| : Split |enterprise_d|, yielding
//                                    |enterprise_d_saucer|.
// * |HistoryTime(2)|               : Remove |enterprise|, merge |enterprise_d|
//                                    and |enterprise_d_saucer|.
// * |HistoryTime(2 + half_a_step)| : Physics bubble ends.
TEST_F(PluginTest, PhysicsBubble) {
  Time const δt = 0.02 * Second;
  Time const a_while = 10 * δt;
  Time const half_a_step = 250 * δt;
  Time const ε_δt = 0.1 * δt;
  EXPECT_THAT(half_a_step, Eq(plugin_->Δt() / 2));
  GUID const enterprise = "NCC-1701";
  GUID const enterprise_d = "NCC-1701-D";
  GUID const enterprise_d_saucer = "NCC-1701-D (saucer)";
  auto const make_enterprise_whole_ship = []() {
    return std::make_pair(
        PartId(0U),
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin,
                Velocity<World>({1 * Metre / Second,
                                 0 * Metre / Second,
                                 0 * Metre / Second})),
            1 * Kilogram,
            Vector<Acceleration, World>({1 * Metre / Pow<2>(Second),
                                         0 * Metre / Pow<2>(Second),
                                         0 * Metre / Pow<2>(Second)})));
  };
  auto const make_enterprise_d_engineering_section = []() {
    return std::make_pair(
        PartId(1U),
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin +
                    Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre}),
                Velocity<World>({0 * Metre / Second,
                                 1 * Metre / Second,
                                 0 * Metre / Second})),
            3 * Kilogram,
            Vector<Acceleration, World>({1 * Metre / Pow<2>(Second),
                                         0 * Metre / Pow<2>(Second),
                                         0 * Metre / Pow<2>(Second)})));
  };
  auto const make_enterprise_d_saucer_section = []() {
    return std::make_pair(
        PartId(2U),
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin -
                    Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre}),
                -Velocity<World>({0 * Metre / Second,
                                  1 * Metre / Second,
                                  0 * Metre / Second})),
            5 * Kilogram,
            Vector<Acceleration, World>({1 * Metre / Pow<2>(Second),
                                         0 * Metre / Pow<2>(Second),
                                         0 * Metre / Pow<2>(Second)})));
  };
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  Angle const planetarium_rotation = 42 * Radian;
  std::size_t expected_number_of_new_off_rails_vessels = 0;
  std::size_t expected_number_of_dirty_old_on_rails_vessels = 0;
  std::size_t expected_number_of_clean_old_vessels = 0;
  std::vector<IdAndOwnedPart> parts;
  bool expect_to_have_physics_bubble = true;
  bool expect_intrinsic_acceleration = false;
  InsertVessel(enterprise, &expected_number_of_new_off_rails_vessels);
  --expected_number_of_new_off_rails_vessels;
  for (int step = 0; step < 10; ++step) {
    for (Instant t = HistoryTime(step) + 2 * δt;
         t <= HistoryTime(step + 1);
         t += δt) {
      // Keep our vessels.  Make sure we're not inserting new ones.
      if (t < HistoryTime(0) + a_while + ε_δt) {
        KeepVessel(enterprise);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_whole_ship()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise, std::move(parts));
      } else if (t > HistoryTime(0) + a_while + ε_δt &&
                 t < HistoryTime(1) + half_a_step + ε_δt) {
        KeepVessel(enterprise_d);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_d_engineering_section()));
        parts.emplace_back(std::move(make_enterprise_d_saucer_section()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
        KeepVessel(enterprise);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_whole_ship()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise, std::move(parts));
      } else if (t > HistoryTime(1) + half_a_step - ε_δt &&
                 t < HistoryTime(2) + ε_δt) {
        KeepVessel(enterprise_d);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_d_engineering_section()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
        KeepVessel(enterprise_d_saucer);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_d_saucer_section()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise_d_saucer,
                                              std::move(parts));
        KeepVessel(enterprise);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_whole_ship()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise, std::move(parts));
      } else if (t > HistoryTime(2) - ε_δt &&
                 t < HistoryTime(2) + half_a_step - ε_δt) {
        KeepVessel(enterprise_d);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_d_engineering_section()));
        parts.emplace_back(std::move(make_enterprise_d_saucer_section()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
      } else {
        KeepVessel(enterprise_d);
      }
      if (AbsoluteError(t - HistoryTime(2), half_a_step) < ε_δt) {
        ++expected_number_of_dirty_old_on_rails_vessels;
        expect_to_have_physics_bubble = false;
      }
      expect_intrinsic_acceleration &= expect_to_have_physics_bubble;
      // Called to compute the prolongations and advance the unsynchronized
      // histories.
      if (expect_intrinsic_acceleration) {
        EXPECT_CALL(
            *n_body_system_,
            Integrate(Ref(plugin_->prolongation_integrator()), t,
                      plugin_->Δt(), 0, true,
                      AllOf(
                          SizeIs(bodies_.size() +
                                 expected_number_of_clean_old_vessels +
                                 expected_number_of_new_off_rails_vessels +
                                 expected_number_of_dirty_old_on_rails_vessels +
                                 (expect_to_have_physics_bubble ? 1 : 0)),
                          Contains(HasNonvanishingIntrinsicAccelerationAt(t)))))
            .RetiresOnSaturation();
      } else {
        EXPECT_CALL(
            *n_body_system_,
            Integrate(Ref(plugin_->prolongation_integrator()), t,
                      plugin_->Δt(), 0, true,
                      SizeIs(bodies_.size() +
                             expected_number_of_clean_old_vessels +
                             expected_number_of_new_off_rails_vessels +
                             expected_number_of_dirty_old_on_rails_vessels +
                             (expect_to_have_physics_bubble ? 1 : 0))))
            .RetiresOnSaturation();
      }
      plugin_->AdvanceTime(t, planetarium_rotation);
      if (expect_to_have_physics_bubble) {
        plugin_->BubbleDisplacementCorrection(World::origin);
        plugin_->BubbleVelocityCorrection(SolarSystem::kSaturn);
      }
      if (AbsoluteError(t - HistoryTime(0), a_while) < ε_δt) {
        InsertVessel(enterprise_d, &expected_number_of_new_off_rails_vessels);
        --expected_number_of_new_off_rails_vessels;
      } else if (AbsoluteError(t - HistoryTime(1), half_a_step) < ε_δt) {
        InsertVessel(enterprise_d_saucer,
                     &expected_number_of_new_off_rails_vessels);
        --expected_number_of_new_off_rails_vessels;
      }
      expect_intrinsic_acceleration = expect_to_have_physics_bubble;
    }
    // Keep the vessels for the history-advancing step.
    if (step <= 0) {
      KeepVessel(enterprise);
      parts.clear();
      parts.emplace_back(std::move(make_enterprise_whole_ship()));
      plugin_->AddVesselToNextPhysicsBubble(enterprise, std::move(parts));
    }
    KeepVessel(enterprise_d);
    if (step <= 1) {
      parts.clear();
      parts.emplace_back(std::move(make_enterprise_d_engineering_section()));
      parts.emplace_back(std::move(make_enterprise_d_saucer_section()));
      plugin_->AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
    }
    // Called to advance the synchronized histories.
    EXPECT_CALL(
        *n_body_system_,
        Integrate(Ref(plugin_->history_integrator()),
                      HistoryTime(step + 1) + δt,
                      plugin_->Δt(), 0, false,
                      SizeIs(bodies_.size() +
                             expected_number_of_clean_old_vessels)))
        .WillOnce(AppendTimeToTrajectories<5>(HistoryTime(step + 1)))
        .RetiresOnSaturation();
    if (expected_number_of_new_off_rails_vessels > 0 ||
        expected_number_of_dirty_old_on_rails_vessels > 0 ||
        expect_to_have_physics_bubble) {
      // Called to synchronize the new histories.
      EXPECT_CALL(
          *n_body_system_,
          Integrate(Ref(plugin_->prolongation_integrator()),
                        HistoryTime(step + 1),
                        plugin_->Δt(), 0, true,
                        SizeIs(bodies_.size() +
                               expected_number_of_new_off_rails_vessels +
                               expected_number_of_dirty_old_on_rails_vessels +
                               (expect_to_have_physics_bubble ? 1 : 0))))
          .WillOnce(AppendTimeToTrajectories<5>(HistoryTime(step + 1)))
          .RetiresOnSaturation();
    }
    expected_number_of_clean_old_vessels +=
        expected_number_of_new_off_rails_vessels +
        expected_number_of_dirty_old_on_rails_vessels;
    expected_number_of_new_off_rails_vessels = 0;
    expected_number_of_dirty_old_on_rails_vessels = 0;
    // Called to compute the prolongations.
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->prolongation_integrator()),
                          HistoryTime(step + 1) + δt, plugin_->Δt(), 0, true,
                          SizeIs(bodies_.size() +
                                     expected_number_of_clean_old_vessels +
                                     (expect_to_have_physics_bubble ? 1 : 0))))
        .RetiresOnSaturation();
    plugin_->AdvanceTime(HistoryTime(step + 1) + δt, planetarium_rotation);
    if (expect_to_have_physics_bubble) {
      plugin_->BubbleDisplacementCorrection(World::origin);
      plugin_->BubbleVelocityCorrection(SolarSystem::kSaturn);
    }
  }
}

TEST_F(PluginTest, UpdateCelestialHierarchy) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    plugin_->UpdateCelestialHierarchy(index, SolarSystem::kSun);
  }
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    RelativeDegreesOfFreedom<ICRFJ2000Ecliptic> const from_parent =
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom() -
        solar_system_->trajectories()[SolarSystem::kSun]->
            last().degrees_of_freedom();
    EXPECT_THAT(
        from_parent,
        Componentwise(
            AlmostEquals(looking_glass_.Inverse()(
                plugin_->CelestialFromParent(index).displacement()), 1, 5056),
            AlmostEquals(looking_glass_.Inverse()(
                plugin_->CelestialFromParent(index).velocity()), 1, 936)));
  }
}

TEST_F(PluginTest, BodyCentredNonrotatingRenderingIntegration) {
  GUID const satellite = "satellite";
  // This is an integration test, so we need a plugin that will actually
  // integrate.
  Plugin plugin(initial_time_,
                SolarSystem::kSun,
                sun_gravitational_parameter_,
                planetarium_rotation_);
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    Index const parent_index = SolarSystem::parent(index);
    RelativeDegreesOfFreedom<AliceSun> const from_parent= looking_glass_(
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom() -
        solar_system_->trajectories()[parent_index]->
            last().degrees_of_freedom());
    plugin.InsertCelestial(index,
                           bodies_[index]->gravitational_parameter(),
                           parent_index,
                           from_parent);
  }
  plugin.EndInitialization();
  plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  plugin.SetVesselStateOffset(satellite,
                              RelativeDegreesOfFreedom<AliceSun>(
                                  satellite_initial_displacement_,
                                  satellite_initial_velocity_));
  not_null<std::unique_ptr<
      Transforms<Barycentric, Rendering, Barycentric>>> const geocentric =
          plugin.NewBodyCentredNonRotatingTransforms(SolarSystem::kEarth);
  // We'll check that our orbit is rendered as circular (actually, we only check
  // that it is rendered within a thin spherical shell around the Earth).
  Length perigee = std::numeric_limits<double>::infinity() * Metre;
  Length apogee = -std::numeric_limits<double>::infinity() * Metre;
  Permutation<AliceSun, World> const alice_sun_to_world =
      Permutation<AliceSun, World>(Permutation<AliceSun, World>::XZY);
  Time const δt_long = 10 * Minute;
#if !defined(_DEBUG)
  Time const δt_short = 0.02 * Second;
  Instant t = initial_time_ + δt_short;
  // Exercise #267 by having small time steps at the beginning of the trajectory
  // that are not synchronized with those of the Earth.
  for (; t < initial_time_ + δt_long; t += δt_short) {
    plugin.AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  }
#else
  Instant t = initial_time_ + δt_long;
  plugin.AdvanceTime(t, 0 * Radian);  // This ensures the history is nonempty.
  plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  t += δt_long;
#endif
  for (; t < initial_time_ + 12 * Hour; t += δt_long) {
    plugin.AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
    // We give the sun an arbitrary nonzero velocity in |World|.
    Position<World> const sun_world_position =
        World::origin + Velocity<World>(
            { 0.1 * AstronomicalUnit / Hour,
             -1.0 * AstronomicalUnit / Hour,
              0.0 * AstronomicalUnit / Hour}) * (t - initial_time_);
    RenderedTrajectory<World> const rendered_trajectory =
        plugin.RenderedVesselTrajectory(satellite,
                                        geocentric.get(),
                                        sun_world_position);
    Position<World> const earth_world_position =
        sun_world_position + alice_sun_to_world(
            plugin.CelestialFromParent(SolarSystem::kEarth).displacement());
    for (auto const segment : rendered_trajectory) {
      Length const l_min =
          std::min((segment.begin - earth_world_position).Norm(),
                   (segment.end - earth_world_position).Norm());
      Length const l_max =
          std::max((segment.begin - earth_world_position).Norm(),
                   (segment.end - earth_world_position).Norm());
      perigee = std::min(perigee, l_min);
      apogee = std::max(apogee, l_max);
    }
    // Check continuity.
    for (std::size_t i = 0; i + 1 < rendered_trajectory.size(); ++i) {
      EXPECT_THAT(rendered_trajectory[i].end,
                  Eq(rendered_trajectory[i + 1].begin));
    }
    EXPECT_THAT(Abs(apogee - perigee), Lt(1.1 * Metre));
  }
}

TEST_F(PluginTest, BarycentricRotatingRenderingIntegration) {
  GUID const satellite = "satellite";
  // This is an integration test, so we need a plugin that will actually
  // integrate.
  Plugin plugin(initial_time_,
                SolarSystem::kSun,
                sun_gravitational_parameter_,
                planetarium_rotation_);
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    Index const parent_index = SolarSystem::parent(index);
    RelativeDegreesOfFreedom<AliceSun> const from_parent = looking_glass_(
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom() -
        solar_system_->trajectories()[parent_index]->
            last().degrees_of_freedom());
    plugin.InsertCelestial(index,
                           bodies_[index]->gravitational_parameter(),
                           parent_index,
                           from_parent);
  }
  plugin.EndInitialization();
  plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  // A vessel at the Lagrange point L₅.
  RelativeDegreesOfFreedom<ICRFJ2000Ecliptic> const from_the_earth_to_the_moon =
      solar_system_->trajectories()[SolarSystem::kMoon]->
          last().degrees_of_freedom() -
      solar_system_->trajectories()[SolarSystem::kEarth]->
          last().degrees_of_freedom();
  Displacement<ICRFJ2000Ecliptic> const from_the_earth_to_l5 =
      from_the_earth_to_the_moon.displacement() / 2 -
          Normalize(from_the_earth_to_the_moon.velocity()) *
              from_the_earth_to_the_moon.displacement().Norm() * Sqrt(3) / 2;
  Velocity<ICRFJ2000Ecliptic> const initial_velocity =
      Rotation<ICRFJ2000Ecliptic, ICRFJ2000Ecliptic>(
          π / 3 * Radian,
          Wedge(from_the_earth_to_the_moon.velocity(),
                from_the_earth_to_the_moon.displacement()))(
              from_the_earth_to_the_moon.velocity());
  plugin.SetVesselStateOffset(satellite,
                              looking_glass_(
                                  RelativeDegreesOfFreedom<ICRFJ2000Ecliptic>(
                                      from_the_earth_to_l5,
                                      initial_velocity)));
  not_null<std::unique_ptr<
      Transforms<Barycentric, Rendering, Barycentric>>> const
      earth_moon_barycentric =
          plugin.NewBarycentricRotatingTransforms(SolarSystem::kEarth,
                                                  SolarSystem::kMoon);
  Permutation<AliceSun, World> const alice_sun_to_world =
      Permutation<AliceSun, World>(Permutation<AliceSun, World>::XZY);
  Time const δt_long = 1 * Hour;
#if defined(_DEBUG)
  Time const duration = 1 * Day;
  Instant t = initial_time_ + δt_long;
#else
  Time const δt_short = 0.02 * Second;
  Time const duration = 20 * Day;
  Instant t = initial_time_ + δt_short;
  // Exercise #267 by having small time steps at the beginning of the trajectory
  // that are not synchronized with those of the Earth and Moon.
  for (; t < initial_time_ + δt_long; t += δt_short) {
    plugin.AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  }
#endif
  for (; t < initial_time_ + duration; t += δt_long) {
    plugin.AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  }
  plugin.AdvanceTime(t,
                     1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
  plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  // We give the sun an arbitrary nonzero velocity in |World|.
  Position<World> const sun_world_position =
      World::origin + Velocity<World>(
          { 0.1 * AstronomicalUnit / Hour,
           -1.0 * AstronomicalUnit / Hour,
            0.0 * AstronomicalUnit / Hour}) * (t - initial_time_);
  RenderedTrajectory<World> const rendered_trajectory =
      plugin.RenderedVesselTrajectory(satellite,
                                      earth_moon_barycentric.get(),
                                      sun_world_position);
  Position<World> const earth_world_position =
      sun_world_position + alice_sun_to_world(
          plugin.CelestialFromParent(SolarSystem::kEarth).displacement());
  Position<World> const moon_world_position =
      earth_world_position + alice_sun_to_world(
          plugin.CelestialFromParent(SolarSystem::kMoon).displacement());
  Length const earth_moon =
      (moon_world_position - earth_world_position).Norm();
  for (auto const segment : rendered_trajectory) {
    Length const satellite_earth =
        (segment.begin - earth_world_position).Norm();
    Length const satellite_moon =
        (segment.begin - moon_world_position).Norm();
    EXPECT_THAT(RelativeError(earth_moon, satellite_earth), Lt(0.0907));
    EXPECT_THAT(RelativeError(earth_moon, satellite_moon), Lt(0.131));
    EXPECT_THAT(RelativeError(satellite_moon, satellite_earth), Lt(0.148));
  }
  // Check continuity.
  for (std::size_t i = 0; i + 1 < rendered_trajectory.size(); ++i) {
    EXPECT_THAT(rendered_trajectory[i].end,
                Eq(rendered_trajectory[i + 1].begin));
  }
#if !defined(_DEBUG)
  // Check that there are no spikes in the rendered trajectory, i.e., that three
  // consecutive points form a sufficiently flat triangle.  This tests issue
  // #256.
  for (std::size_t i = 0; i + 2 < rendered_trajectory.size(); ++i) {
    EXPECT_THAT(
        (rendered_trajectory[i].begin - rendered_trajectory[i + 1].end).Norm(),
        Gt(((rendered_trajectory[i].begin -
                 rendered_trajectory[i + 1].begin).Norm() +
            (rendered_trajectory[i].end -
                 rendered_trajectory[i + 1].end).Norm()) / 1.5)) << i;
  }
#endif
}

// Checks that we correctly predict a full circular orbit around a massive body
// with unit gravitational parameter at unit distance, in 8 steps.  Since
// predictions are only computed on |AdvanceTime()|, we advance time by a small
// amount.
TEST_F(PluginTest, Prediction) {
  GUID const satellite = "satellite";
  Index const celestial = 0;
  int const n = 8;
  Plugin plugin(Instant(),
                celestial,
                SIUnit<GravitationalParameter>(),
                0 * Radian);
  plugin.EndInitialization();
  EXPECT_TRUE(plugin.InsertOrKeepVessel(satellite, celestial));
  auto transforms = plugin.NewBodyCentredNonRotatingTransforms(celestial);
  plugin.SetVesselStateOffset(
      satellite,
      {Displacement<AliceSun>({1 * Metre, 0 * Metre, 0 * Metre}),
       Velocity<AliceSun>(
           {0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second})});
  plugin.set_predicted_vessel(satellite);
  plugin.set_prediction_length(2 * π * Second);
  plugin.set_prediction_step(2 * π / n * Second);
  plugin.AdvanceTime(Instant(1e-10 * Second), 0 * Radian);
  RenderedTrajectory<World> rendered_prediction =
      plugin.RenderedPrediction(transforms.get(), World::origin);
  EXPECT_EQ(n, rendered_prediction.size());
  Angle const α = 2 * π * Radian / n;
  for (int k = 0; k < n; ++k) {
    EXPECT_THAT(
        RelativeError(rendered_prediction[k].begin - World::origin,
                      Displacement<World>({Cos(k * α) * Metre,
                                           0 * Metre,
                                           Sin(k * α) * Metre})),
        Lt(0.011));
    EXPECT_THAT(
        RelativeError(rendered_prediction[k].end - World::origin,
                      Displacement<World>({Cos((k + 1) * α) * Metre,
                                           0 * Metre,
                                           Sin((k + 1) * α) * Metre})),
        Lt(0.011));
  }
  plugin.clear_predicted_vessel();
}

TEST_F(PluginTest, NavBall) {
  // Create a plugin with planetarium rotation 0.
  auto plugin = Plugin(initial_time_,
                       SolarSystem::kSun,
                       sun_gravitational_parameter_,
                       0 * Radian);
  not_null<std::unique_ptr<
      Transforms<Barycentric, Rendering, Barycentric>>> const heliocentric =
          plugin.NewBodyCentredNonRotatingTransforms(SolarSystem::kSun);
  Vector<double, World> x({1, 0, 0});
  Vector<double, World> y({0, 1, 0});
  Vector<double, World> z({0, 0, 1});
  auto nav_ball = plugin.NavBall(heliocentric.get(), World::origin);
  EXPECT_THAT(AbsoluteError(-z, nav_ball(World::origin)(x)),
              Lt(2 * std::numeric_limits<double>::epsilon()));
  EXPECT_THAT(AbsoluteError(y, nav_ball(World::origin)(y)),
              Lt(std::numeric_limits<double>::epsilon()));
  EXPECT_THAT(AbsoluteError(x, nav_ball(World::origin)(z)),
              Lt(2 * std::numeric_limits<double>::epsilon()));
}

}  // namespace ksp_plugin
}  // namespace principia
