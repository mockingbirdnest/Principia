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
#include "physics/mock_ephemeris.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system.hpp"

namespace principia {

using geometry::Bivector;
using geometry::Permutation;
using geometry::Trivector;
using physics::MockEphemeris;
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
// This parameter must be a |Trajectories|, |time| must be an |Instant|.
ACTION_TEMPLATE(AppendTimeToTrajectories,
                HAS_1_TEMPLATE_PARAMS(int, k),
                AND_1_VALUE_PARAMS(time)) {
  for (auto trajectory :
          static_cast<std::vector<not_null<Trajectory<Barycentric>*>>>(
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
  std::map<Index, ContinuousTrajectory<Barycentric>> trajectories_;
  std::unique_ptr<StrictMock<MockEphemeris<Barycentric>>> mock_ephemeris_;

 public:
  TestablePlugin(Instant const& initial_time,
                 Index const sun_index,
                 GravitationalParameter const& sun_gravitational_parameter,
                 Angle const& planetarium_rotation)
      : Plugin(initial_time,
               planetarium_rotation),
        mock_ephemeris_(
            std::make_unique<StrictMock<MockEphemeris<Barycentric>>>()) {}

  // We need to override |EndInitialization in order to create a |MockEphemeris|
  // rather than an |Ephemeris|, and in order to fill in the continuous
  // trajectories ourselves.
  void EndInitialization() override {
    initializing_.Flop();
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
    std::vector<DegreesOfFreedom<Barycentric>> initial_state;
    auto bodies_it = bodies_->begin();
    for (auto const& index_state : *initial_state_) {
      EXPECT_EQ(index_state.first, bodies_it->first);
      auto const inserted =
          trajectories_.emplace(std::piecewise_construct,
                                std::forward_as_tuple(index_state.first),
                                std::forward_as_tuple(45 * Minute,
                                                      1 * Milli(Metre)));
      EXPECT_TRUE(inserted.second);
      for (int i = 0; i < 9; ++i) {
        inserted.first->second.Append(
            current_time_ + i * 45 * Minute,
            {index_state.second.position() +
                 i * 45 * Minute * index_state.second.velocity(),
             index_state.second.velocity()});
      }
      ++bodies_it;
    }
    bodies_.reset();
    initial_state_.reset();
    ephemeris_ = std::move(mock_ephemeris_);
    for (auto const& index_celestial : celestials_) {
      auto const& index = index_celestial.first;
      auto& celestial = *index_celestial.second;
      celestial.set_trajectory(&FindOrDie(trajectories_, index));
    }
  }

  Time const& Δt() const {
    return Δt_;
  }

  StrictMock<MockEphemeris<Barycentric>>* mock_ephemeris() const {
    return mock_ephemeris_.get();
  }
};

class PluginTest : public testing::Test {
 protected:
  PluginTest()
      : looking_glass_(Permutation<ICRFJ2000Ecliptic, AliceSun>::XZY),
        solar_system_(SolarSystem::AtСпутник1Launch(
            SolarSystem::Accuracy::kMajorBodiesOnly)),
        bodies_(solar_system_->massive_bodies()),
        initial_time_(42 * Second),
        sun_gravitational_parameter_(
            bodies_[SolarSystem::kSun]->gravitational_parameter()),
        planetarium_rotation_(1 * Radian),
        plugin_(make_not_null_unique<TestablePlugin>(
                    initial_time_,
                    SolarSystem::kSun,
                    sun_gravitational_parameter_,
                    planetarium_rotation_)) {
    mock_ephemeris_ = plugin_->mock_ephemeris();
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
  plugin_->InsertSun(SolarSystem::kSun, sun_gravitational_parameter_);
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
  Instant HistoryTime(Instant const sync_time, int const step) {
    return sync_time + step * plugin_->Δt();
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
                    not_null<std::size_t*> const number_of_new_vessels,
                    Instant const& time) {
    bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                      SolarSystem::kEarth);
    EXPECT_TRUE(inserted) << guid;
    EXPECT_CALL(*mock_ephemeris_, Prolong(time)).RetiresOnSaturation();
    plugin_->SetVesselStateOffset(guid,
                                  RelativeDegreesOfFreedom<AliceSun>(
                                      satellite_initial_displacement_,
                                      satellite_initial_velocity_));
    ++*number_of_new_vessels;
  }

  Permutation<ICRFJ2000Ecliptic, AliceSun> looking_glass_;
  StrictMock<MockEphemeris<Barycentric>>* mock_ephemeris_;
  not_null<std::unique_ptr<SolarSystem>> solar_system_;
  SolarSystem::Bodies bodies_;
  Instant initial_time_;
  GravitationalParameter sun_gravitational_parameter_;
  Angle planetarium_rotation_;

  not_null<std::unique_ptr<TestablePlugin>> plugin_;

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
                    planetarium_rotation_);
  plugin->InsertSun(SolarSystem::kSun, sun_gravitational_parameter_);
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    Index const parent_index = SolarSystem::parent(index);
    RelativeDegreesOfFreedom<AliceSun> const from_parent = looking_glass_(
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

  Instant const& sync_time = initial_time_ + 1 * Second;
  // Sync.
  plugin->AdvanceTime(sync_time, Angle());

  // Add a handful of points to the history and then forget some of them.  This
  // is the most convenient way to check that forgetting works as expected.
  plugin->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  plugin->AdvanceTime(HistoryTime(sync_time, 3), Angle());
  plugin->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  plugin->AdvanceTime(HistoryTime(sync_time, 6), Angle());
  plugin->ForgetAllHistoriesBefore(HistoryTime(sync_time, 3));

  serialization::Plugin message;
  plugin->WriteToMessage(&message);
  plugin = Plugin::ReadFromMessage(message);
  serialization::Plugin second_message;
  plugin->WriteToMessage(&second_message);
  EXPECT_EQ(message.SerializeAsString(), second_message.SerializeAsString());
  EXPECT_EQ(bodies_.size(), message.celestial_size());

  EXPECT_FALSE(message.celestial(0).has_parent_index());
  EXPECT_EQ(message.celestial(0).index(), message.celestial(1).parent_index());

  EXPECT_EQ(
      HistoryTime(sync_time, 3),
      Instant::ReadFromMessage(message.ephemeris().trajectory(0).first_time()));

  EXPECT_EQ(1, message.vessel_size());
  EXPECT_EQ(SolarSystem::kEarth, message.vessel(0).parent_index());
  EXPECT_TRUE(message.vessel(0).vessel().has_history_and_prolongation());
  auto const& vessel_0_history =
      message.vessel(0).vessel().history_and_prolongation().history();
  EXPECT_EQ(3, vessel_0_history.timeline_size());
  EXPECT_EQ((HistoryTime(sync_time, 4) - Instant()) / (1 * Second),
            vessel_0_history.timeline(0).instant().scalar().magnitude());
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
                        1, 278784),
                    AlmostEquals(looking_glass_.Inverse()(
                            plugin_->CelestialFromParent(index).velocity()),
                        0, 1643885)));
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

TEST_F(PluginDeathTest, SunError) {
  EXPECT_DEATH({
    plugin_->InsertSun(42, bodies_.front()->gravitational_parameter());
    plugin_->InsertSun(43, bodies_.front()->gravitational_parameter());
  }, "sun_ == nullptr");
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
    EXPECT_CALL(*mock_ephemeris_, Prolong(initial_time_));
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


TEST_F(PluginDeathTest, ForgetAllHistoriesBeforeError) {
  EXPECT_DEATH({
    Instant const t = initial_time_ + 100 * Second;
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->AdvanceTime(t, Angle());
    plugin_->ForgetAllHistoriesBefore(t);
  }, "Check failed: t < history_time_");
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
  EXPECT_CALL(*mock_ephemeris_, Prolong(initial_time_));
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
// This now checks that we do nothing but prolong, since there are no vessels.
TEST_F(PluginTest, AdvanceTimeWithCelestialsOnly) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  Time const δt = 0.02 * Second;
  Angle const planetarium_rotation = 42 * Radian;
  for (int step = 0; step < 10; ++step) {
    for (Instant t = HistoryTime(initial_time_, step) + 2 * δt;
         t <= HistoryTime(initial_time_, step + 1);
         t += δt) {
      EXPECT_CALL(*mock_ephemeris_, Prolong(t)).RetiresOnSaturation();
      plugin_->AdvanceTime(t, planetarium_rotation);
    }
    EXPECT_CALL(*mock_ephemeris_,
                Prolong(HistoryTime(initial_time_, step + 1) + δt))
        .RetiresOnSaturation();
    plugin_->AdvanceTime(HistoryTime(initial_time_, step + 1) + δt,
                         planetarium_rotation);
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
  Time const δt = 1 * Second;
  Time const a_while = 2 * δt;
  Time const half_a_step = 5 * δt;
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
  InSequence in_sequence;
  EXPECT_FALSE(plugin_->has_vessel(enterprise));
  Instant sync_time = initial_time_;
  InsertVessel(enterprise,
               &expected_number_of_new_vessels,
               HistoryTime(sync_time, 0));
  EXPECT_TRUE(plugin_->has_vessel(enterprise));
  EXPECT_CALL(*mock_ephemeris_,
              Prolong(HistoryTime(sync_time, 0) + δt))
      .RetiresOnSaturation();
  // Called to compute the prolongation.
  EXPECT_CALL(*mock_ephemeris_,
              FlowWithAdaptiveStep(_, _, _, _, HistoryTime(sync_time, 0) + δt))
      .RetiresOnSaturation();
  // Therer are no old vessels, this will make all the new vessels old
  // instantly without performing any syncing computations.
  plugin_->AdvanceTime(HistoryTime(sync_time, 0) + δt, planetarium_rotation);
  expected_number_of_old_vessels = expected_number_of_new_vessels;
  expected_number_of_new_vessels = 0;
  sync_time = HistoryTime(sync_time, 0) + δt;
  for (int step = 0; step < 10; ++step) {
    for (Instant t = HistoryTime(sync_time, step) + 2 * δt;
         t < HistoryTime(sync_time, step + 1);
         t += δt) {
      // Keep our vessels.  Make sure we're not inserting new ones.
      if (step <= 3) {
        KeepVessel(enterprise);
      }
      if (t > HistoryTime(sync_time, 0) + a_while + ε_δt) {
        KeepVessel(enterprise_d);
      }
      if (t > HistoryTime(sync_time, 1) + a_while + ε_δt) {
        KeepVessel(stargazer);
      }
      if (t > HistoryTime(sync_time, 1) + a_while + ε_δt &&
          t < HistoryTime(sync_time, 1) + half_a_step - ε_δt) {
        KeepVessel(bradbury);
      } else if (AbsoluteError(t - HistoryTime(sync_time, 1),
                               half_a_step) < ε_δt) {
        // We will be removing |bradbury| in this step.
        --expected_number_of_new_vessels;
      }
      if (step > 2) {
        KeepVessel(constantinople);
      }
      EXPECT_CALL(*mock_ephemeris_, Prolong(t)).RetiresOnSaturation();
      for (int i = 0;
           i < expected_number_of_new_vessels + expected_number_of_old_vessels;
           ++i) {
        // Called to compute the prolongations and advance the unsynchronized
        // histories.
        EXPECT_CALL(*mock_ephemeris_,
                    FlowWithAdaptiveStep(_, _, _, _, t)).RetiresOnSaturation();
      }
      plugin_->AdvanceTime(t, planetarium_rotation);
      if (AbsoluteError(t - HistoryTime(sync_time, 0), a_while) < ε_δt) {
        InsertVessel(enterprise_d, &expected_number_of_new_vessels, t);
      } else if (AbsoluteError(t - HistoryTime(sync_time, 1), a_while) < ε_δt) {
        InsertVessel(stargazer, &expected_number_of_new_vessels, t);
        InsertVessel(bradbury, &expected_number_of_new_vessels, t);
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
    EXPECT_CALL(*mock_ephemeris_,
                Prolong(HistoryTime(sync_time, step + 1) + δt))
        .RetiresOnSaturation();
    if (expected_number_of_old_vessels == 0) {
      expected_number_of_old_vessels = expected_number_of_new_vessels;
      expected_number_of_new_vessels = 0;
    } else {
      // Called to advance the synchronized histories.
      EXPECT_CALL(*mock_ephemeris_,
                  FlowWithFixedStep(SizeIs(expected_number_of_old_vessels),
                                    10 * Second,
                                    HistoryTime(sync_time, step + 1) + δt))
          .WillOnce(
               AppendTimeToTrajectories<0>(HistoryTime(sync_time, step + 1)))
          .RetiresOnSaturation();
    }
    for (int i = 0; i < expected_number_of_new_vessels; ++i) {
      // Called to synchronize the new histories.
      EXPECT_CALL(*mock_ephemeris_,
                  FlowWithAdaptiveStep(_, _, _, _,
                                       HistoryTime(sync_time, step + 1)))
          .RetiresOnSaturation();
    }
    expected_number_of_old_vessels += expected_number_of_new_vessels;
    expected_number_of_new_vessels = 0;
    for (int i = 0; i < expected_number_of_old_vessels; ++i) {
      // Called to compute the prolongations.
      EXPECT_CALL(*mock_ephemeris_,
                  FlowWithAdaptiveStep(_, _, _, _,
                                       HistoryTime(sync_time, step + 1) + δt))
          .RetiresOnSaturation();
    }
    plugin_->AdvanceTime(HistoryTime(sync_time, step + 1) + δt,
                         planetarium_rotation);
    if (step == 2) {
      InsertVessel(constantinople,
                   &expected_number_of_new_vessels,
                   HistoryTime(sync_time, step + 1) + δt);
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
  Time const δt = 1 * Second;
  Time const a_while = 2 * δt;
  Time const half_a_step = 5 * δt;
  Time const ε_δt = 0.1 * δt;
  EXPECT_THAT(half_a_step, Eq(plugin_->Δt() / 2));
  GUID const enterprise = "NCC-1701";
  GUID const enterprise_d = "NCC-1701-D";
  GUID const enterprise_d_saucer = "NCC-1701-D (saucer)";
  GUID const constantinople = "NCC-43622";
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
  Instant sync_time = initial_time_;
  // The synchronizer vessel.
  InsertVessel(constantinople,
               &expected_number_of_clean_old_vessels,
               HistoryTime(sync_time, 0));
  // Called to compute the prolongation.
  EXPECT_CALL(*mock_ephemeris_,
              FlowWithAdaptiveStep(_, _, _, _, HistoryTime(sync_time, 0) + δt))
      .RetiresOnSaturation();
  // There are no old vessels, this will make all the new vessels old instantly
  // without performing any syncing computations.
  EXPECT_CALL(*mock_ephemeris_,
              Prolong(HistoryTime(sync_time, 0) + δt))
      .RetiresOnSaturation();
  plugin_->AdvanceTime(HistoryTime(sync_time, 0) + δt,
                       planetarium_rotation);
  sync_time = HistoryTime(sync_time, 0) + δt;

  InsertVessel(enterprise,
               &expected_number_of_new_off_rails_vessels,
               HistoryTime(sync_time, 0));
  --expected_number_of_new_off_rails_vessels;
  for (int step = 0; step < 10; ++step) {
    for (Instant t = HistoryTime(sync_time, step) + 2 * δt;
         t < HistoryTime(sync_time, step + 1);
         t += δt) {
      // Keep the synchronizer.
      KeepVessel(constantinople);
      // Keep our vessels.  Make sure we're not inserting new ones.
      if (t < HistoryTime(sync_time, 0) + a_while + ε_δt) {
        KeepVessel(enterprise);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_whole_ship()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise, std::move(parts));
      } else if (t > HistoryTime(sync_time, 0) + a_while + ε_δt &&
                 t < HistoryTime(sync_time, 1) + half_a_step + ε_δt) {
        KeepVessel(enterprise_d);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_d_engineering_section()));
        parts.emplace_back(std::move(make_enterprise_d_saucer_section()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
        KeepVessel(enterprise);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_whole_ship()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise, std::move(parts));
      } else if (t > HistoryTime(sync_time, 1) + half_a_step - ε_δt &&
                 t < HistoryTime(sync_time, 2) + ε_δt) {
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
      } else if (t > HistoryTime(sync_time, 2) - ε_δt &&
                 t < HistoryTime(sync_time, 2) + half_a_step - ε_δt) {
        KeepVessel(enterprise_d);
        parts.clear();
        parts.emplace_back(std::move(make_enterprise_d_engineering_section()));
        parts.emplace_back(std::move(make_enterprise_d_saucer_section()));
        plugin_->AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
      } else {
        KeepVessel(enterprise_d);
      }
      if (AbsoluteError(t - HistoryTime(sync_time, 2), half_a_step) < ε_δt) {
        ++expected_number_of_dirty_old_on_rails_vessels;
        expect_to_have_physics_bubble = false;
      }
      expect_intrinsic_acceleration &= expect_to_have_physics_bubble;
      // Called to compute the prolongations and advance the unsynchronized
      // histories.
      for (int i = 0;
           i < expected_number_of_clean_old_vessels +
               expected_number_of_new_off_rails_vessels +
               expected_number_of_dirty_old_on_rails_vessels +
               (expect_to_have_physics_bubble ? 1 : 0);
           ++i) {
      EXPECT_CALL(*mock_ephemeris_,
                  FlowWithAdaptiveStep(_, _, _, _, t))
            .RetiresOnSaturation();
      }
      EXPECT_CALL(*mock_ephemeris_, Prolong(t)).RetiresOnSaturation();
      plugin_->AdvanceTime(t, planetarium_rotation);
      if (expect_to_have_physics_bubble) {
        plugin_->BubbleDisplacementCorrection(World::origin);
        plugin_->BubbleVelocityCorrection(SolarSystem::kSaturn);
      }
      if (AbsoluteError(t - HistoryTime(sync_time, 0), a_while) < ε_δt) {
        InsertVessel(enterprise_d,
                     &expected_number_of_new_off_rails_vessels,
                     t);
        --expected_number_of_new_off_rails_vessels;
      } else if (AbsoluteError(t - HistoryTime(sync_time, 1),
                               half_a_step) < ε_δt) {
        InsertVessel(enterprise_d_saucer,
                     &expected_number_of_new_off_rails_vessels,
                     t);
        --expected_number_of_new_off_rails_vessels;
      }
      expect_intrinsic_acceleration = expect_to_have_physics_bubble;
    }
    // Keep the vessels for the history-advancing step.
    // Keep the synchronizer.
    KeepVessel(constantinople);
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
        *mock_ephemeris_,
        FlowWithFixedStep(SizeIs(expected_number_of_clean_old_vessels),
                          10 * Second,
                          HistoryTime(sync_time, step + 1) + δt))
        .WillOnce(
             AppendTimeToTrajectories<0>(HistoryTime(sync_time, step + 1)))
        .RetiresOnSaturation();
     for (int i = 0;
          i < expected_number_of_new_off_rails_vessels +
              expected_number_of_dirty_old_on_rails_vessels +
              (expect_to_have_physics_bubble ? 1 : 0);
          ++i) {
      // Called to synchronize the new histories.
      EXPECT_CALL(
          *mock_ephemeris_,
          FlowWithAdaptiveStep(_, _, _, _, HistoryTime(sync_time, step + 1)))
          .RetiresOnSaturation();
    }
    expected_number_of_clean_old_vessels +=
        expected_number_of_new_off_rails_vessels +
        expected_number_of_dirty_old_on_rails_vessels;
    expected_number_of_new_off_rails_vessels = 0;
    expected_number_of_dirty_old_on_rails_vessels = 0;
    // Called to compute the prolongations.
    for (int i = 0;
         i < expected_number_of_clean_old_vessels +
             (expect_to_have_physics_bubble ? 1 : 0);
         ++i) {
      // Called to compute the prolongations.
      EXPECT_CALL(
          *mock_ephemeris_,
          FlowWithAdaptiveStep(_, _, _, _,
                               HistoryTime(sync_time, step + 1) + δt))
          .RetiresOnSaturation();
    }
    EXPECT_CALL(*mock_ephemeris_,
                Prolong(HistoryTime(sync_time, step + 1) + δt))
        .RetiresOnSaturation();
    plugin_->AdvanceTime(HistoryTime(sync_time, step + 1) + δt,
                         planetarium_rotation);
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
                plugin_->CelestialFromParent(index).velocity()), 1, 1643885)));
  }
}

TEST_F(PluginTest, BodyCentredNonrotatingRenderingIntegration) {
  GUID const satellite = "satellite";
  // This is an integration test, so we need a plugin that will actually
  // integrate.
  Plugin plugin(initial_time_,
                planetarium_rotation_);
  plugin.InsertSun(SolarSystem::kSun, sun_gravitational_parameter_);
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
  plugin.SetVesselStateOffset(satellite,
                              RelativeDegreesOfFreedom<AliceSun>(
                                  satellite_initial_displacement_,
                                  satellite_initial_velocity_));
  not_null<std::unique_ptr<RenderingTransforms>> const geocentric =
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
                planetarium_rotation_);
  plugin.InsertSun(SolarSystem::kSun, sun_gravitational_parameter_);
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
  not_null<std::unique_ptr<RenderingTransforms>> const earth_moon_barycentric =
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
// with unit gravitational parameter at unit distance.  Since predictions are
// only computed on |AdvanceTime()|, we advance time by a small amount.
TEST_F(PluginTest, Prediction) {
  GUID const satellite = "satellite";
  Index const celestial = 0;
  Plugin plugin(Instant(),
                0 * Radian);
  plugin.InsertSun(celestial, SIUnit<GravitationalParameter>());
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
  plugin.set_prediction_length_tolerance(1 * Milli(Metre));
  plugin.set_prediction_speed_tolerance(1 * Milli(Metre) / Second);
  plugin.AdvanceTime(Instant(1e-10 * Second), 0 * Radian);
  RenderedTrajectory<World> rendered_prediction =
      plugin.RenderedPrediction(transforms.get(), World::origin);
  EXPECT_EQ(14, rendered_prediction.size());
  for (int i = 0; i < rendered_prediction.size(); ++i) {
    auto const& segment = rendered_prediction[i];
    EXPECT_THAT(
        AbsoluteError((segment.begin - World::origin).Norm(), 1 * Metre),
        Lt(0.5 * Milli(Metre)));
    EXPECT_THAT(AbsoluteError((segment.end - World::origin).Norm(), 1 * Metre),
                Lt(0.5 * Milli(Metre)));
    if (i >= 5) {
      EXPECT_THAT(
          AbsoluteError((segment.begin - World::origin).Norm(), 1 * Metre),
          Gt(0.1 * Milli(Metre)));
      EXPECT_THAT(
          AbsoluteError((segment.end - World::origin).Norm(), 1 * Metre),
          Gt(0.1 * Milli(Metre)));
    }
  }
  EXPECT_THAT(
      AbsoluteError(rendered_prediction.back().end - World::origin,
                    Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre})),
      AllOf(Gt(2 * Milli(Metre)), Lt(3 * Milli(Metre))));
  plugin.clear_predicted_vessel();
}

TEST_F(PluginTest, Navball) {
  // Create a plugin with planetarium rotation 0.
  Plugin plugin(initial_time_,
                0 * Radian);
  plugin.InsertSun(SolarSystem::kSun, sun_gravitational_parameter_);
  plugin.EndInitialization();
  not_null<std::unique_ptr<RenderingTransforms>> const heliocentric =
          plugin.NewBodyCentredNonRotatingTransforms(SolarSystem::kSun);
  Vector<double, World> x({1, 0, 0});
  Vector<double, World> y({0, 1, 0});
  Vector<double, World> z({0, 0, 1});
  auto navball = plugin.Navball(heliocentric.get(), World::origin);
  EXPECT_THAT(AbsoluteError(-z, navball(World::origin)(x)),
              Lt(2 * std::numeric_limits<double>::epsilon()));
  EXPECT_THAT(AbsoluteError(y, navball(World::origin)(y)),
              Lt(std::numeric_limits<double>::epsilon()));
  EXPECT_THAT(AbsoluteError(x, navball(World::origin)(z)),
              Lt(2 * std::numeric_limits<double>::epsilon()));
}

TEST_F(PluginTest, Frenet) {
  // Create a plugin with planetarium rotation 0.
  Plugin plugin(initial_time_,
                0 * Radian);
  plugin.InsertSun(SolarSystem::kEarth,
                   bodies_[SolarSystem::kEarth]->gravitational_parameter());
  plugin.EndInitialization();
  Permutation<AliceSun, World> const alice_sun_to_world =
      Permutation<AliceSun, World>(Permutation<AliceSun, World>::XZY);
  GUID const satellite = "satellite";
  plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  plugin.SetVesselStateOffset(satellite,
                              RelativeDegreesOfFreedom<AliceSun>(
                                  satellite_initial_displacement_,
                                  satellite_initial_velocity_));
  Vector<double, World> t = alice_sun_to_world(
                                Normalize(satellite_initial_velocity_));
  not_null<std::unique_ptr<RenderingTransforms>> const geocentric =
          plugin.NewBodyCentredNonRotatingTransforms(SolarSystem::kEarth);
  EXPECT_THAT(plugin.VesselTangent(satellite, geocentric.get()),
              AlmostEquals(t, 7));
}

#if 0
TEST_F(PluginTest, SerializationCompatibility) {
  serialization::Multivector message;

  Vector<Length, Barycentric> const v({-1 * Metre, 2 * Metre, 3 * Metre});
  v.WriteToMessage(&message);
  message.mutable_frame()->set_tag(serialization::Frame::OLD_BARYCENTRIC);
  Vector<Length, Barycentric> const w =
      Vector<Length, Barycentric>::ReadFromMessage(message);
  Vector<Length, Barycentric> const expected_w(
      {-1 * Metre, 3 * Metre, 2 * Metre});
  EXPECT_EQ(expected_w, w);

  Bivector<Length, Barycentric> const b({4 * Metre, 5 * Metre, -6 * Metre});
  b.WriteToMessage(&message);
  message.mutable_frame()->set_tag(serialization::Frame::OLD_BARYCENTRIC);
  Bivector<Length, Barycentric> const c =
      Bivector<Length, Barycentric>::ReadFromMessage(message);
  Bivector<Length, Barycentric> const expected_c(
      {-4 * Metre, 6 * Metre, -5 * Metre});
  EXPECT_EQ(expected_c, c);

  Trivector<Length, Barycentric> const t(-7 * Metre);
  t.WriteToMessage(&message);
  message.mutable_frame()->set_tag(serialization::Frame::OLD_BARYCENTRIC);
  Trivector<Length, Barycentric> const u =
      Trivector<Length, Barycentric>::ReadFromMessage(message);
  Trivector<Length, Barycentric> const expected_u(7 * Metre);
  EXPECT_EQ(expected_u, u);
}
#endif

}  // namespace ksp_plugin
}  // namespace principia
