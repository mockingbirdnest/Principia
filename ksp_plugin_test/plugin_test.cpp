
#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/macros.hpp"
#include "geometry/identity.hpp"
#include "geometry/permutation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/kepler_orbit.hpp"
#include "physics/mock_dynamic_frame.hpp"
#include "physics/mock_ephemeris.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::Bivector;
using geometry::Permutation;
using geometry::Trivector;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MockDynamicFrame;
using physics::MockEphemeris;
using physics::RigidTransformation;
using quantities::Abs;
using quantities::ArcTan;
using quantities::Cos;
using quantities::Sin;
using quantities::Sqrt;
using quantities::si::Day;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Minute;
using quantities::si::Newton;
using quantities::si::Radian;
using quantities::si::AstronomicalUnit;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystemFactory;
using ::testing::AllOf;
using ::testing::AnyNumber;
using ::testing::Contains;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::InSequence;
using ::testing::Le;
using ::testing::Lt;
using ::testing::Ref;
using ::testing::Return;
using ::testing::SizeIs;
using ::testing::StrictMock;
using ::testing::_;

namespace ksp_plugin {

namespace {

int const kNotABody = 1729;

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

ACTION(AppendToDiscreteTrajectories) {
  for (auto const& trajectory : arg0) {
    trajectory->Append(arg3, {Barycentric::origin, Velocity<Barycentric>()});
  }
}

ACTION(AppendToDiscreteTrajectory) {
  arg0->Append(arg5, {Barycentric::origin, Velocity<Barycentric>()});
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

  Time const& Δt() const {
    return Δt_;
  }

  StrictMock<MockEphemeris<Barycentric>>* mock_ephemeris() const {
    return mock_ephemeris_.get();
  }

  Rotation<AliceSun, Barycentric> InversePlanetariumRotation() {
    return PlanetariumRotation().Inverse();
  }

 private:
  // We override this part of initialization in order to create a
  // |MockEphemeris| rather than an |Ephemeris|, and in order to fill in the
  // continuous trajectories ourselves.
  void InitializeEphemerisAndSetCelestialTrajectories() override {
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
    std::vector<DegreesOfFreedom<Barycentric>> initial_state;
    auto bodies_it = absolute_initialization_->bodies.begin();
    for (auto const& index_state : absolute_initialization_->initial_state) {
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
    absolute_initialization_ = std::experimental::nullopt;
    ephemeris_ = std::move(mock_ephemeris_);
    for (auto const& index_celestial : celestials_) {
      auto const& index = index_celestial.first;
      auto& celestial = *index_celestial.second;
      celestial.set_trajectory(&FindOrDie(trajectories_, index));
    }
  }
};

class PluginTest : public testing::Test {
 protected:
  PluginTest()
      : solar_system_(SolarSystemFactory::AtСпутник1Launch(
            SolarSystemFactory::Accuracy::kMajorBodiesOnly)),
        initial_time_(Instant() + 42 * Second),
        sun_gravitational_parameter_(
            solar_system_->gravitational_parameter(
                SolarSystemFactory::name(SolarSystemFactory::kSun))),
        sun_mean_radius_(
            solar_system_->mean_radius(
                SolarSystemFactory::name(SolarSystemFactory::kSun))),
        planetarium_rotation_(1 * Radian),
        plugin_(make_not_null_unique<TestablePlugin>(
                    initial_time_,
                    SolarSystemFactory::kSun,
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
        Sqrt(solar_system_->gravitational_parameter(
                 SolarSystemFactory::name(SolarSystemFactory::kEarth)) /
                 satellite_initial_displacement_.Norm()) * unit_tangent;
  }

  void InsertAllSolarSystemBodies() {
    for (int index = SolarSystemFactory::kSun;
         index <= SolarSystemFactory::kLastMajorBody;
         ++index) {
      std::experimental::optional<Index> parent_index;
      if (index != SolarSystemFactory::kSun) {
        parent_index = SolarSystemFactory::parent(index);
      }
      std::string const name = SolarSystemFactory::name(index);
      plugin_->InsertCelestialAbsoluteCartesian(
          index,
          parent_index,
          id_icrf_barycentric_(
              solar_system_->initial_state(SolarSystemFactory::name(index))),
          make_not_null_unique<MassiveBody>(
              solar_system_->gravitational_parameter(name)));
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
    bool const inserted = plugin_->InsertOrKeepVessel(
                              guid, SolarSystemFactory::kEarth);
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
    bool const inserted = plugin_->InsertOrKeepVessel(
                              guid, SolarSystemFactory::kEarth);
    EXPECT_TRUE(inserted) << guid;
    EXPECT_CALL(*mock_ephemeris_, Prolong(time)).RetiresOnSaturation();
    plugin_->SetVesselStateOffset(guid,
                                  RelativeDegreesOfFreedom<AliceSun>(
                                      satellite_initial_displacement_,
                                      satellite_initial_velocity_));
    ++*number_of_new_vessels;
  }

  static RigidMotion<ICRFJ2000Equator, Barycentric> const id_icrf_barycentric_;
  StrictMock<MockEphemeris<Barycentric>>* mock_ephemeris_;
  not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>> solar_system_;
  Instant const initial_time_;
  GravitationalParameter const sun_gravitational_parameter_;
  Length const sun_mean_radius_;
  Angle planetarium_rotation_;

  not_null<std::unique_ptr<TestablePlugin>> plugin_;

  // These initial conditions will yield a low circular orbit around Earth.
  Displacement<AliceSun> satellite_initial_displacement_;
  Velocity<AliceSun> satellite_initial_velocity_;
};

RigidMotion<ICRFJ2000Equator, Barycentric> const
    PluginTest::id_icrf_barycentric_(
        RigidTransformation<ICRFJ2000Equator, Barycentric>(
            ICRFJ2000Equator::origin,
            Barycentric::origin,
            OrthogonalMap<ICRFJ2000Equator, Barycentric>::Identity()),
        AngularVelocity<ICRFJ2000Equator>(),
        Velocity<ICRFJ2000Equator>());

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
  plugin->InsertSun(SolarSystemFactory::kSun,
                    sun_gravitational_parameter_,
                    sun_mean_radius_);
  for (int index = SolarSystemFactory::kSun + 1;
       index <= SolarSystemFactory::kLastMajorBody;
       ++index) {
    std::string const name = SolarSystemFactory::name(index);
    Index const parent_index = SolarSystemFactory::parent(index);
    std::string const parent_name = SolarSystemFactory::name(parent_index);
    RelativeDegreesOfFreedom<Barycentric> const state_vectors =
        Identity<ICRFJ2000Equator, Barycentric>()(
            solar_system_->initial_state(name) -
            solar_system_->initial_state(parent_name));
    Instant const t;
    auto body = make_not_null_unique<MassiveBody>(
        solar_system_->gravitational_parameter(name));
    KeplerianElements<Barycentric> elements = KeplerOrbit<Barycentric>(
        /*primary=*/MassiveBody(
            solar_system_->gravitational_parameter(parent_name)),
        /*secondary=*/*body,
        state_vectors,
        /*epoch=*/t).elements_at_epoch();
    elements.semimajor_axis = std::experimental::nullopt;
    plugin->InsertCelestialJacobiKeplerian(index,
                                           parent_index,
                                           elements,
                                           std::move(body));
  }
  plugin->EndInitialization();
  plugin->InsertOrKeepVessel(satellite, SolarSystemFactory::kEarth);
  plugin->SetVesselStateOffset(satellite,
                               RelativeDegreesOfFreedom<AliceSun>(
                                   satellite_initial_displacement_,
                                   satellite_initial_velocity_));

  Instant const& sync_time = initial_time_ + 1 * Second;
  // Sync.
  plugin->AdvanceTime(sync_time, Angle());

  // Add a handful of points to the history and then forget some of them.  This
  // is the most convenient way to check that forgetting works as expected.
  plugin->InsertOrKeepVessel(satellite, SolarSystemFactory::kEarth);
  plugin->AdvanceTime(HistoryTime(sync_time, 3), Angle());
  plugin->InsertOrKeepVessel(satellite, SolarSystemFactory::kEarth);
  plugin->AdvanceTime(HistoryTime(sync_time, 6), Angle());
  plugin->ForgetAllHistoriesBefore(HistoryTime(sync_time, 3));

  plugin->CreateFlightPlan(satellite, HistoryTime(sync_time, 7), 4 * Kilogram);
  plugin->SetPlottingFrame(plugin->NewBodyCentredNonRotatingNavigationFrame(
      SolarSystemFactory::kSun + 1));

  serialization::Plugin message;
  plugin->WriteToMessage(&message);
  plugin = Plugin::ReadFromMessage(message);
  serialization::Plugin second_message;
  plugin->WriteToMessage(&second_message);
  EXPECT_EQ(message.SerializeAsString(), second_message.SerializeAsString())
      << "FIRST\n" << message.DebugString()
      << "SECOND\n" << second_message.DebugString();;
  EXPECT_EQ(SolarSystemFactory::kLastMajorBody - SolarSystemFactory::kSun + 1,
            message.celestial_size());

  EXPECT_FALSE(message.celestial(0).has_parent_index());
  EXPECT_EQ(message.celestial(0).index(), message.celestial(1).parent_index());

  EXPECT_EQ(
      HistoryTime(sync_time, 3),
      Instant::ReadFromMessage(message.ephemeris().trajectory(0).first_time()));

  EXPECT_EQ(1, message.vessel_size());
  EXPECT_EQ(SolarSystemFactory::kEarth, message.vessel(0).parent_index());
  EXPECT_TRUE(message.vessel(0).vessel().has_flight_plan());
  EXPECT_TRUE(message.vessel(0).vessel().has_history_and_prolongation());
  auto const& vessel_0_history =
      message.vessel(0).vessel().history_and_prolongation().history();
#if defined(WE_LOVE_228)
  EXPECT_EQ(2, vessel_0_history.timeline_size());
  EXPECT_EQ((HistoryTime(sync_time, 3) - Instant()) / (1 * Second),
            vessel_0_history.timeline(0).instant().scalar().magnitude());
  EXPECT_EQ((HistoryTime(sync_time, 6) - Instant()) / (1 * Second),
            vessel_0_history.timeline(1).instant().scalar().magnitude());
#else
  EXPECT_EQ(3, vessel_0_history.timeline_size());
  EXPECT_EQ((HistoryTime(sync_time, 4) - Instant()) / (1 * Second),
            vessel_0_history.timeline(0).instant().scalar().magnitude());
#endif
  EXPECT_FALSE(message.bubble().has_current());
  EXPECT_TRUE(message.has_plotting_frame());
  EXPECT_TRUE(message.plotting_frame().HasExtension(
      serialization::BodyCentredNonRotatingDynamicFrame::
          body_centred_non_rotating_dynamic_frame));
  EXPECT_EQ(SolarSystemFactory::kSun + 1,
            message.plotting_frame().GetExtension(
                serialization::BodyCentredNonRotatingDynamicFrame::
                    body_centred_non_rotating_dynamic_frame).centre());
}

TEST_F(PluginTest, Initialization) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  EXPECT_CALL(*mock_ephemeris_, Prolong(_)).Times(AnyNumber());
  for (int index = SolarSystemFactory::kSun + 1;
       index <= SolarSystemFactory::kLastMajorBody;
       ++index) {
    auto const to_icrf = id_icrf_barycentric_.orthogonal_map().Inverse() *
                         plugin_->InversePlanetariumRotation().Forget();
    Index const parent_index = SolarSystemFactory::parent(index);
    RelativeDegreesOfFreedom<ICRFJ2000Equator> const from_parent =
        solar_system_->initial_state(SolarSystemFactory::name(index)) -
        solar_system_->initial_state(SolarSystemFactory::name(parent_index));
    EXPECT_THAT(from_parent,
                Componentwise(
                    AlmostEquals(to_icrf(plugin_->CelestialFromParent(index)
                                             .displacement()),
                                 0, 42380),
                    AlmostEquals(
                        to_icrf(plugin_->CelestialFromParent(index).velocity()),
                        74, 1475468)))
        << SolarSystemFactory::name(index);
  }
}

TEST_F(PluginTest, HierarchicalInitialization) {
  // e, i, Ω, ω, and mean anomaly are 0.
  KeplerianElements<Barycentric> elements;

  // We construct a system as follows, inserting the bodies in the order
  // S0, P1, P2, M3.
  // |<1 m>|     |<1 m>|
  // 2     1     1     2
  //   |<   7/3 m   >|
  // S0    P2    M3    P1
  plugin_->InsertSun(0, 2 * SIUnit<GravitationalParameter>(), 1 * Metre);
  elements.semimajor_axis = 7.0 / 3.0 * Metre;
  plugin_->InsertCelestialJacobiKeplerian(
      /*celestial_index=*/1,
      /*parent_index=*/0,
      elements,
      make_not_null_unique<MassiveBody>(2 * SIUnit<GravitationalParameter>()));
  elements.semimajor_axis = 1 * Metre;
  plugin_->InsertCelestialJacobiKeplerian(
      /*celestial_index=*/2,
      /*parent_index=*/0,
      elements,
      make_not_null_unique<MassiveBody>(1 * SIUnit<GravitationalParameter>()));
  elements.mean_anomaly = π * Radian;
  plugin_->InsertCelestialJacobiKeplerian(
      /*celestial_index=*/3,
      /*parent_index=*/1,
      elements,
      make_not_null_unique<MassiveBody>(1 * SIUnit<GravitationalParameter>()));
  plugin_->EndInitialization();
  EXPECT_CALL(*mock_ephemeris_, Prolong(_)).Times(AnyNumber());
  EXPECT_THAT(plugin_->CelestialFromParent(1).displacement().Norm(),
              AlmostEquals(3.0 * Metre, 1));
  EXPECT_THAT(plugin_->CelestialFromParent(2).displacement().Norm(),
              Eq(1 * Metre));
  EXPECT_THAT(plugin_->CelestialFromParent(3).displacement().Norm(),
              Eq(1 * Metre));
}

TEST_F(PluginDeathTest, SunError) {
  EXPECT_DEATH({
    plugin_->InsertSun(42, sun_gravitational_parameter_, sun_mean_radius_);
    plugin_->InsertSun(43, sun_gravitational_parameter_, sun_mean_radius_);
  }, "!hierarchical_initialization_");
}

TEST_F(PluginDeathTest, UpdateCelestialHierarchyError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->UpdateCelestialHierarchy(SolarSystemFactory::kSun,
                                      SolarSystemFactory::kPluto);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->UpdateCelestialHierarchy(kNotABody, SolarSystemFactory::kPluto);
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->UpdateCelestialHierarchy(SolarSystemFactory::kSun, kNotABody);
  }, "Map key not found");
}

TEST_F(PluginDeathTest, InsertOrKeepVesselError) {
  GUID const guid = "Syrio Forel";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::kSun);
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
    plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::kSun);
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

TEST_F(PluginTest, ForgetAllHistoriesBeforeWithFlightPlan) {
  GUID const guid = "Test Satellite";
  Instant const t = initial_time_ + 100 * Second;

  auto* const mock_dynamic_frame =
      new MockDynamicFrame<Barycentric, Navigation>();
  EXPECT_CALL(*mock_ephemeris_, Prolong(_)).Times(AnyNumber());
  EXPECT_CALL(*mock_ephemeris_, FlowWithAdaptiveStep(_, _, _, _, _, _))
      .WillRepeatedly(AppendToDiscreteTrajectory());
  EXPECT_CALL(*mock_ephemeris_, FlowWithFixedStep(_, _, _, _))
      .WillRepeatedly(AppendToDiscreteTrajectories());
  EXPECT_CALL(*mock_ephemeris_, ForgetBefore(_)).Times(1);
  EXPECT_CALL(*mock_dynamic_frame, ToThisFrameAtTime(_))
      .WillRepeatedly(Return(
          RigidMotion<Barycentric, Navigation>(
              RigidTransformation<Barycentric, Navigation>::Identity(),
              AngularVelocity<Barycentric>(),
              Velocity<Barycentric>())));
  EXPECT_CALL(*mock_dynamic_frame, FrenetFrame(_, _))
      .WillRepeatedly(Return(
          MockDynamicFrame<Barycentric, Navigation>::Rot::Identity()));

  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();

  plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::kEarth);
  plugin_->SetVesselStateOffset(guid,
                                RelativeDegreesOfFreedom<AliceSun>(
                                    satellite_initial_displacement_,
                                    satellite_initial_velocity_));
  auto const satellite = plugin_->GetVessel(guid);

  Instant const& sync_time = initial_time_ + 1 * Second;
  plugin_->AdvanceTime(sync_time, Angle());
  plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::kEarth);
  plugin_->AdvanceTime(HistoryTime(sync_time, 3), Angle());

  auto const burn = [this, mock_dynamic_frame, sync_time]() -> Burn {
    return {/*thrust=*/1 * Newton,
            /*specific_impulse=*/1 * Newton * Second / Kilogram,
            std::unique_ptr<MockDynamicFrame<Barycentric, Navigation>>(
                mock_dynamic_frame),
            /*initial_time=*/HistoryTime(sync_time, 4),
            Velocity<Frenet<Navigation>>(
                {1 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second})};
  };
  plugin_->CreateFlightPlan(guid,
                            /*final_time=*/HistoryTime(sync_time, 8),
                            /*initial_mass=*/1 * Kilogram);
  satellite->flight_plan()->Append(burn());

  plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::kEarth);
  plugin_->AdvanceTime(HistoryTime(sync_time, 6), Angle());
  plugin_->ForgetAllHistoriesBefore(HistoryTime(sync_time, 1));
  EXPECT_LE(satellite->history().Begin().time(),
            satellite->flight_plan()->initial_time());
  EXPECT_EQ(1 * Newton, satellite->flight_plan()->GetManœuvre(0).thrust());
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
    plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::kSun);
    plugin_->VesselFromParent(guid);
  }, "not given an initial state");
}

TEST_F(PluginDeathTest, CelestialFromParentError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->CelestialFromParent(SolarSystemFactory::kEarth);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->CelestialFromParent(kNotABody);
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->CelestialFromParent(SolarSystemFactory::kSun);
  }, "is the sun");
}

TEST_F(PluginTest, VesselInsertionAtInitialization) {
  GUID const guid = "Test Satellite";
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                    SolarSystemFactory::kEarth);
  EXPECT_TRUE(inserted);
  EXPECT_CALL(*mock_ephemeris_, Prolong(initial_time_)).Times(AnyNumber());
  plugin_->SetVesselStateOffset(guid,
                                RelativeDegreesOfFreedom<AliceSun>(
                                    satellite_initial_displacement_,
                                    satellite_initial_velocity_));
  EXPECT_THAT(plugin_->VesselFromParent(guid),
              Componentwise(
                  AlmostEquals(satellite_initial_displacement_, 14496),
                  AlmostEquals(satellite_initial_velocity_, 3)));
}

TEST_F(PluginTest, UpdateCelestialHierarchy) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  EXPECT_CALL(*mock_ephemeris_, Prolong(_)).Times(AnyNumber());
  for (int index = SolarSystemFactory::kSun + 1;
       index <= SolarSystemFactory::kLastMajorBody;
       ++index) {
    plugin_->UpdateCelestialHierarchy(index, SolarSystemFactory::kSun);
  }
  for (int index = SolarSystemFactory::kSun + 1;
       index <= SolarSystemFactory::kLastMajorBody;
       ++index) {
    auto const to_icrf = id_icrf_barycentric_.orthogonal_map().Inverse() *
                         plugin_->InversePlanetariumRotation().Forget();
    RelativeDegreesOfFreedom<ICRFJ2000Equator> const from_parent =
        solar_system_->initial_state(SolarSystemFactory::name(index)) -
        solar_system_->initial_state(
            SolarSystemFactory::name(SolarSystemFactory::kSun));
    // All these worlds are fine -- except Triton.
    // Attempt no computation there.
    if (index == SolarSystemFactory::kTriton) {
      EXPECT_THAT(
          from_parent,
          Componentwise(
              AlmostEquals(
                  to_icrf(plugin_->CelestialFromParent(index).displacement()),
                  2),
              AlmostEquals(
                  to_icrf(plugin_->CelestialFromParent(index).velocity()),
                  20155840)))
          << SolarSystemFactory::name(index);
    } else {
      EXPECT_THAT(
          from_parent,
          Componentwise(
              AlmostEquals(
                  to_icrf(plugin_->CelestialFromParent(index).displacement()),
                  0, 13),
              AlmostEquals(
                  to_icrf(plugin_->CelestialFromParent(index).velocity()),
                  74, 1475468)))
          << SolarSystemFactory::name(index);
    }
  }
}
TEST_F(PluginTest, Navball) {
  // Create a plugin with planetarium rotation 0.
  Plugin plugin(initial_time_,
                0 * Radian);
  plugin.InsertSun(SolarSystemFactory::kSun,
                   sun_gravitational_parameter_,
                   sun_mean_radius_);
  plugin.EndInitialization();
  not_null<std::unique_ptr<NavigationFrame>> navigation_frame =
      plugin.NewBodyCentredNonRotatingNavigationFrame(SolarSystemFactory::kSun);
  not_null<const NavigationFrame*> const navigation_frame_copy =
      navigation_frame.get();
  plugin.SetPlottingFrame(std::move(navigation_frame));
  EXPECT_EQ(navigation_frame_copy, plugin.GetPlottingFrame());
  Vector<double, World> x({1, 0, 0});
  Vector<double, World> y({0, 1, 0});
  Vector<double, World> z({0, 0, 1});
  auto navball = plugin.Navball(World::origin);
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
  plugin.InsertSun(SolarSystemFactory::kEarth,
                   solar_system_->gravitational_parameter(
                       SolarSystemFactory::name(SolarSystemFactory::kEarth)),
                   solar_system_->mean_radius(
                       SolarSystemFactory::name(SolarSystemFactory::kEarth)));
  plugin.EndInitialization();
  Permutation<AliceSun, World> const alice_sun_to_world =
      Permutation<AliceSun, World>(Permutation<AliceSun, World>::XZY);
  GUID const satellite = "satellite";
  plugin.InsertOrKeepVessel(satellite, SolarSystemFactory::kEarth);
  plugin.SetVesselStateOffset(satellite,
                              RelativeDegreesOfFreedom<AliceSun>(
                                  satellite_initial_displacement_,
                                  satellite_initial_velocity_));
  Vector<double, World> t = alice_sun_to_world(
                                Normalize(satellite_initial_velocity_));
  Vector<double, World> n = alice_sun_to_world(
                                Normalize(-satellite_initial_displacement_));
  // World is left-handed, but the Frenet trihedron is right-handed.
  Vector<double, World> b(-geometry::Cross(t.coordinates(), n.coordinates()));
  not_null<std::unique_ptr<NavigationFrame>> const geocentric =
      plugin.NewBodyCentredNonRotatingNavigationFrame(
          SolarSystemFactory::kEarth);
  EXPECT_THAT(plugin.VesselTangent(satellite), AlmostEquals(t, 2));
  EXPECT_THAT(plugin.VesselNormal(satellite), AlmostEquals(n, 3));
  EXPECT_THAT(plugin.VesselBinormal(satellite), AlmostEquals(b, 4));
}

}  // namespace ksp_plugin
}  // namespace principia
