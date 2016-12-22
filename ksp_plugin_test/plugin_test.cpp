
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
#include "base/not_null.hpp"
#include "geometry/identity.hpp"
#include "geometry/permutation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/mock_dynamic_frame.hpp"
#include "physics/mock_ephemeris.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system_factory.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

using astronomy::ICRFJ2000Equator;
using base::FindOrDie;
using base::make_not_null_unique;
using base::not_null;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::Identity;
using geometry::Permutation;
using geometry::Trivector;
using integrators::McLachlanAtela1992Order5Optimal;
using physics::ContinuousTrajectory;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MockDynamicFrame;
using physics::MockEphemeris;
using physics::RigidMotion;
using physics::RigidTransformation;
using physics::SolarSystem;
using quantities::Abs;
using quantities::AngularFrequency;
using quantities::ArcTan;
using quantities::Cos;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Sin;
using quantities::SIUnit;
using quantities::Sqrt;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Minute;
using quantities::si::Newton;
using quantities::si::Radian;
using quantities::si::AstronomicalUnit;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystemFactory;
using testing_utilities::VanishesBefore;
using ::testing::AllOf;
using ::testing::AnyNumber;
using ::testing::Contains;
using ::testing::DoAll;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::InSequence;
using ::testing::Le;
using ::testing::Lt;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::SetArgPointee;
using ::testing::SizeIs;
using ::testing::StrictMock;
using ::testing::_;

namespace {

int const not_a_body = 1729;

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
    trajectory->Append(arg2, {Barycentric::origin, Velocity<Barycentric>()});
  }
}

ACTION(AppendToDiscreteTrajectory) {
  arg0->Append(arg2, {Barycentric::origin, Velocity<Barycentric>()});
}

}  // namespace

class TestablePlugin : public Plugin {
 public:
  TestablePlugin(Instant const& game_epoch,
                 Instant const& solar_system_epoch,
                 Angle const& planetarium_rotation)
      : Plugin(game_epoch, solar_system_epoch, planetarium_rotation),
        // The |mock_ephemeris_| has to be created early so that we can write
        // expectations before |EndInitialization| has been called.
        owned_mock_ephemeris_(
            std::make_unique<MockEphemeris<Barycentric>>()),
        mock_ephemeris_(owned_mock_ephemeris_.get()) {}

  Time const& Δt() const {
    return history_parameters_.step();
  }

  MockEphemeris<Barycentric>& mock_ephemeris() const {
    return *mock_ephemeris_;
  }

  Rotation<AliceSun, Barycentric> InversePlanetariumRotation() {
    return PlanetariumRotation().Inverse();
  }

  not_null<ContinuousTrajectory<Barycentric> const*> trajectory(
      Index const index) const {
    return trajectories_.at(index).get();
  }

 protected:
  // We override this part of initialization in order to create a
  // |MockEphemeris| rather than an |Ephemeris|.
  std::unique_ptr<Ephemeris<Barycentric>> NewEphemeris(
      std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies,
      std::vector<DegreesOfFreedom<Barycentric>> const& initial_state,
      Instant const& initial_time,
      Length const& fitting_tolerance,
      Ephemeris<Barycentric>::FixedStepParameters const& parameters) override {
    // Own the bodies.
    bodies_ = std::move(bodies);

    // Construct continuous trajectories for the bodies with some data
    Index index = 0;
    for (auto const& degrees_of_freedom : initial_state) {
      trajectories_.push_back(
          base::make_not_null_unique<ContinuousTrajectory<Barycentric>>(
              45 * Minute, 1 * Milli(Metre)));
      auto const& trajectory = trajectories_.back();
      for (int i = 0; i < 9; ++i) {
        trajectory->Append(
            current_time_ + i * 45 * Minute,
            {degrees_of_freedom.position() +
                 i * 45 * Minute * degrees_of_freedom.velocity(),
             degrees_of_freedom.velocity()});
      }

      // Make sure that the |trajectory| member does the right thing.
      ON_CALL(*mock_ephemeris_, trajectory(bodies_[index].get()))
          .WillByDefault(Return(trajectory.get()));

      ++index;
    }

    // Return the mock ephemeris.  We squirelled away a pointer in
    // |mock_ephemeris_|.
    return std::move(owned_mock_ephemeris_);
  }

 private:
  // Someone has to own the bodies; normally it would be the ephemeris, but
  // since we mock it, we keep them here instead.  Similarly for the
  // trajectories.
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies_;
  std::vector<not_null<std::unique_ptr<ContinuousTrajectory<Barycentric>>>>
      trajectories_;
  std::unique_ptr<MockEphemeris<Barycentric>> owned_mock_ephemeris_;
  MockEphemeris<Barycentric>* mock_ephemeris_;
};

class PluginTest : public testing::Test {
 protected:
  PluginTest()
      : body_rotation_(/*mean_radius=*/1 * Metre,
                       /*reference_angle=*/0 * Degree,
                       /*reference_instant=*/astronomy::J2000,
                       /*angular_frequency=*/1 * Radian / Second,
                       /*right_ascension_of_pole=*/0 * Degree,
                       /*declination_of_pole=*/90 * Degree),
        solar_system_(SolarSystemFactory::AtСпутник1Launch(
            SolarSystemFactory::Accuracy::MajorBodiesOnly)),
        initial_time_(Instant() + 42 * Second),
        sun_body_(make_not_null_unique<RotatingBody<Barycentric>>(
            MassiveBody::Parameters(solar_system_->gravitational_parameter(
                SolarSystemFactory::name(SolarSystemFactory::Sun))),
            body_rotation_)),
        planetarium_rotation_(1 * Radian),
        plugin_(make_not_null_unique<TestablePlugin>(
                    initial_time_,
                    initial_time_,
                    planetarium_rotation_)) {
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
                 SolarSystemFactory::name(SolarSystemFactory::Earth)) /
                 satellite_initial_displacement_.Norm()) * unit_tangent;

    // Fill required fields.
    Length{}.WriteToMessage(
        valid_ephemeris_message_.mutable_fitting_tolerance());
    numerics::DoublePrecision<Instant>{}.WriteToMessage(
        valid_ephemeris_message_.mutable_last_state()->mutable_time());
  }

  void InsertAllSolarSystemBodies() {
    for (int index = SolarSystemFactory::Sun;
         index <= SolarSystemFactory::LastMajorBody;
         ++index) {
      std::experimental::optional<Index> parent_index;
      if (index != SolarSystemFactory::Sun) {
        parent_index = SolarSystemFactory::parent(index);
      }
      std::string const name = SolarSystemFactory::name(index);
      plugin_->InsertCelestialAbsoluteCartesian(
          index,
          parent_index,
          id_icrf_barycentric_(
              solar_system_->initial_state(SolarSystemFactory::name(index))),
          SolarSystem<Barycentric>::MakeMassiveBody(
              solar_system_->gravity_model_message(name)));
    }
  }

  // The time of the |step|th history step of |plugin_|.  |HistoryTime(0)| is
  // |initial_time_|.
  Instant HistoryTime(Instant const time, int const step) {
    return time + step * plugin_->Δt();
  }

  // Keeps the vessel with the given |guid| during the next call to
  // |AdvanceTime|.  The vessel must be present.
  void KeepVessel(GUID const& guid) {
    bool const inserted = plugin_->InsertOrKeepVessel(
                              guid, SolarSystemFactory::Earth);
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
                              guid, SolarSystemFactory::Earth);
    EXPECT_TRUE(inserted) << guid;
    EXPECT_CALL(plugin_->mock_ephemeris(), Prolong(time)).RetiresOnSaturation();
    plugin_->SetVesselStateOffset(guid,
                                  RelativeDegreesOfFreedom<AliceSun>(
                                      satellite_initial_displacement_,
                                      satellite_initial_velocity_));
    ++*number_of_new_vessels;
  }

  RotatingBody<Barycentric>::Parameters body_rotation_;
  static RigidMotion<ICRFJ2000Equator, Barycentric> const id_icrf_barycentric_;
  not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>> solar_system_;
  Instant const initial_time_;
  not_null<std::unique_ptr<MassiveBody>> sun_body_;
  Angle planetarium_rotation_;

  not_null<std::unique_ptr<TestablePlugin>> plugin_;

  // These initial conditions will yield a low circular orbit around Earth.
  Displacement<AliceSun> satellite_initial_displacement_;
  Velocity<AliceSun> satellite_initial_velocity_;
  serialization::Ephemeris valid_ephemeris_message_;
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
                    initial_time_,
                    planetarium_rotation_);
  plugin->InsertCelestialJacobiKeplerian(
      SolarSystemFactory::Sun,
      /*parent_index=*/std::experimental::nullopt,
      /*keplerian_elements=*/std::experimental::nullopt,
      std::move(sun_body_));
  for (int index = SolarSystemFactory::Sun + 1;
       index <= SolarSystemFactory::LastMajorBody;
       ++index) {
    std::string const name = SolarSystemFactory::name(index);
    Index const parent_index = SolarSystemFactory::parent(index);
    std::string const parent_name = SolarSystemFactory::name(parent_index);
    RelativeDegreesOfFreedom<Barycentric> const state_vectors =
        Identity<ICRFJ2000Equator, Barycentric>()(
            solar_system_->initial_state(name) -
            solar_system_->initial_state(parent_name));
    Instant const t;
    auto body = make_not_null_unique<RotatingBody<Barycentric>>(
        solar_system_->gravitational_parameter(name),
        body_rotation_);
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
  plugin->InsertOrKeepVessel(satellite, SolarSystemFactory::Earth);
  plugin->SetVesselStateOffset(satellite,
                               RelativeDegreesOfFreedom<AliceSun>(
                                   satellite_initial_displacement_,
                                   satellite_initial_velocity_));

  Time const shift = 1 * Second;
  Instant const time = initial_time_ + shift;
  plugin->AdvanceTime(time, Angle());

  // Add a handful of points to the history and then forget some of them.  This
  // is the most convenient way to check that forgetting works as expected.
  plugin->InsertOrKeepVessel(satellite, SolarSystemFactory::Earth);
  plugin->AdvanceTime(HistoryTime(time, 3), Angle());
  plugin->InsertOrKeepVessel(satellite, SolarSystemFactory::Earth);
  plugin->AdvanceTime(HistoryTime(time, 6), Angle());
  plugin->UpdatePrediction(satellite);
  plugin->ForgetAllHistoriesBefore(HistoryTime(time, 2));

  plugin->CreateFlightPlan(satellite, HistoryTime(time, 7), 4 * Kilogram);
  plugin->SetPlottingFrame(plugin->NewBodyCentredNonRotatingNavigationFrame(
      SolarSystemFactory::Sun + 1));

  serialization::Plugin message;
  plugin->WriteToMessage(&message);
  plugin = Plugin::ReadFromMessage(message);
  serialization::Plugin second_message;
  plugin->WriteToMessage(&second_message);
  EXPECT_EQ(message.SerializeAsString(), second_message.SerializeAsString())
      << "FIRST\n" << message.DebugString()
      << "SECOND\n" << second_message.DebugString();
  EXPECT_EQ(SolarSystemFactory::LastMajorBody - SolarSystemFactory::Sun + 1,
            message.celestial_size());

  EXPECT_FALSE(message.celestial(0).has_parent_index());
  EXPECT_EQ(message.celestial(0).index(), message.celestial(1).parent_index());

  EXPECT_EQ(
      HistoryTime(time, 2),
      Instant::ReadFromMessage(message.ephemeris().trajectory(0).first_time()));

  EXPECT_EQ(1, message.vessel_size());
  EXPECT_EQ(SolarSystemFactory::Earth, message.vessel(0).parent_index());
  EXPECT_TRUE(message.vessel(0).vessel().has_flight_plan());
  EXPECT_TRUE(message.vessel(0).vessel().has_history());
  auto const& vessel_0_history = message.vessel(0).vessel().history();
#if defined(WE_LOVE_228)
  EXPECT_EQ(2, vessel_0_history.timeline_size());
  EXPECT_EQ((HistoryTime(time, 3) - shift - Instant()) / (1 * Second),
            vessel_0_history.timeline(0).instant().scalar().magnitude());
  EXPECT_EQ((HistoryTime(time, 6) - shift - Instant()) / (1 * Second),
            vessel_0_history.timeline(1).instant().scalar().magnitude());
#else
  EXPECT_EQ(3, vessel_0_history.timeline_size());
  EXPECT_EQ((HistoryTime(time, 4) - Instant()) / (1 * Second),
            vessel_0_history.timeline(0).instant().scalar().magnitude());
#endif
  EXPECT_FALSE(message.bubble().has_current());
  EXPECT_TRUE(message.has_plotting_frame());
  EXPECT_TRUE(message.plotting_frame().HasExtension(
      serialization::BodyCentredNonRotatingDynamicFrame::extension));
  EXPECT_EQ(SolarSystemFactory::Sun + 1,
            message.plotting_frame().GetExtension(
                serialization::BodyCentredNonRotatingDynamicFrame::extension).
                    centre());
}

TEST_F(PluginTest, Initialization) {
  InsertAllSolarSystemBodies();
  EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
      .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
  plugin_->EndInitialization();
  EXPECT_CALL(plugin_->mock_ephemeris(), Prolong(_)).Times(AnyNumber());
  for (int index = SolarSystemFactory::Sun + 1;
       index <= SolarSystemFactory::LastMajorBody;
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
  auto sun_body = make_not_null_unique<RotatingBody<Barycentric>>(
      MassiveBody::Parameters(2 * SIUnit<GravitationalParameter>()),
      body_rotation_);
  plugin_->InsertCelestialJacobiKeplerian(
      0,
      /*parent_index=*/std::experimental::nullopt,
      /*keplerian_elements=*/std::experimental::nullopt,
      std::move(sun_body));
  elements.semimajor_axis = 7.0 / 3.0 * Metre;
  plugin_->InsertCelestialJacobiKeplerian(
      /*celestial_index=*/1,
      /*parent_index=*/0,
      elements,
      make_not_null_unique<RotatingBody<Barycentric>>(
          2 * SIUnit<GravitationalParameter>(),
          body_rotation_));
  elements.semimajor_axis = 1 * Metre;
  plugin_->InsertCelestialJacobiKeplerian(
      /*celestial_index=*/2,
      /*parent_index=*/0,
      elements,
      make_not_null_unique<RotatingBody<Barycentric>>(
          1 * SIUnit<GravitationalParameter>(),
          body_rotation_));
  elements.mean_anomaly = π * Radian;
  plugin_->InsertCelestialJacobiKeplerian(
      /*celestial_index=*/3,
      /*parent_index=*/1,
      elements,
      make_not_null_unique<RotatingBody<Barycentric>>(
          1 * SIUnit<GravitationalParameter>(),
          body_rotation_));
  EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
      .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
  plugin_->EndInitialization();
  EXPECT_CALL(plugin_->mock_ephemeris(), Prolong(_)).Times(AnyNumber());
  EXPECT_THAT(plugin_->CelestialFromParent(1).displacement().Norm(),
              AlmostEquals(3.0 * Metre, 3));
  EXPECT_THAT(plugin_->CelestialFromParent(2).displacement().Norm(),
              AlmostEquals(1 * Metre, 2, 3));
  EXPECT_THAT(plugin_->CelestialFromParent(3).displacement().Norm(),
              AlmostEquals(1 * Metre, 2, 3));
}

TEST_F(PluginDeathTest, SunError) {
  EXPECT_DEATH({
      plugin_->InsertCelestialJacobiKeplerian(
          42,
          /*parent_index=*/std::experimental::nullopt,
          /*keplerian_elements=*/std::experimental::nullopt,
          std::move(sun_body_));
      plugin_->InsertCelestialJacobiKeplerian(
          43,
          /*parent_index=*/std::experimental::nullopt,
          /*keplerian_elements=*/std::experimental::nullopt,
          std::move(sun_body_));
  }, ".bool.parent_index == .bool.hierarchical_initialization");
  EXPECT_DEATH({
    KeplerianElements<Barycentric> sun_keplerian_elements;
    sun_keplerian_elements.mean_motion = AngularFrequency();
    plugin_->InsertCelestialJacobiKeplerian(
        43,
        /*parent_index=*/std::experimental::nullopt,
        sun_keplerian_elements,
        std::move(sun_body_));
  }, ".bool.parent_index == .bool.keplerian_elements");
}

TEST_F(PluginDeathTest, UpdateCelestialHierarchyError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->UpdateCelestialHierarchy(SolarSystemFactory::Sun,
                                      SolarSystemFactory::Pluto);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
        .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
    plugin_->EndInitialization();
    plugin_->UpdateCelestialHierarchy(not_a_body, SolarSystemFactory::Pluto);
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
        .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
    plugin_->EndInitialization();
    plugin_->UpdateCelestialHierarchy(SolarSystemFactory::Sun, not_a_body);
  }, "Map key not found");
}

TEST_F(PluginDeathTest, InsertOrKeepVesselError) {
  GUID const guid = "Syrio Forel";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::Sun);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
        .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, not_a_body);
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
    EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
        .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
    plugin_->EndInitialization();
    plugin_->SetVesselStateOffset(guid,
                                  RelativeDegreesOfFreedom<AliceSun>(
                                      satellite_initial_displacement_,
                                      satellite_initial_velocity_));
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
        .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::Sun);
    EXPECT_CALL(plugin_->mock_ephemeris(), Prolong(initial_time_));
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

TEST_F(PluginTest, ForgetAllHistoriesBeforeWithFlightPlan) {
  GUID const guid = "Test Satellite";

  auto* const mock_dynamic_frame =
      new MockDynamicFrame<Barycentric, Navigation>();
  EXPECT_CALL(plugin_->mock_ephemeris(), t_max())
      .WillRepeatedly(Return(Instant()));
  EXPECT_CALL(plugin_->mock_ephemeris(), empty()).WillRepeatedly(Return(false));
  EXPECT_CALL(plugin_->mock_ephemeris(), Prolong(_)).Times(AnyNumber());
  EXPECT_CALL(plugin_->mock_ephemeris(), FlowWithAdaptiveStep(_, _, _, _, _))
      .WillRepeatedly(DoAll(AppendToDiscreteTrajectory(), Return(true)));
  EXPECT_CALL(plugin_->mock_ephemeris(), FlowWithFixedStep(_, _, _, _))
      .WillRepeatedly(AppendToDiscreteTrajectories());
  EXPECT_CALL(plugin_->mock_ephemeris(), planetary_integrator())
      .WillRepeatedly(
          ReturnRef(McLachlanAtela1992Order5Optimal<Position<Barycentric>>()));
  EXPECT_CALL(plugin_->mock_ephemeris(), ForgetBefore(_)).Times(2);
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
  EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
      .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
  plugin_->EndInitialization();

  plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::Earth);
  plugin_->SetVesselStateOffset(guid,
                                RelativeDegreesOfFreedom<AliceSun>(
                                    satellite_initial_displacement_,
                                    satellite_initial_velocity_));
  auto const satellite = plugin_->GetVessel(guid);

  Instant const& time = initial_time_ + 1 * Second;
  plugin_->AdvanceTime(time, Angle());
  plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::Earth);
  plugin_->AdvanceTime(HistoryTime(time, 3), Angle());

  auto const burn = [this, mock_dynamic_frame, time]() -> Burn {
    return {/*thrust=*/1 * Newton,
            /*specific_impulse=*/1 * Newton * Second / Kilogram,
            std::unique_ptr<MockDynamicFrame<Barycentric, Navigation>>(
                mock_dynamic_frame),
            /*initial_time=*/HistoryTime(time, 4),
            Velocity<Frenet<Navigation>>(
                {1 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second})};
  };
  plugin_->CreateFlightPlan(guid,
                            /*final_time=*/HistoryTime(time, 8),
                            /*initial_mass=*/1 * Kilogram);
  satellite->flight_plan().Append(burn());

  plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::Earth);
  plugin_->AdvanceTime(HistoryTime(time, 6), Angle());
  plugin_->ForgetAllHistoriesBefore(HistoryTime(time, 3));
  EXPECT_LE(HistoryTime(time, 3), satellite->flight_plan().initial_time());
  EXPECT_LE(HistoryTime(time, 3), satellite->history().Begin().time());
  EXPECT_EQ(1, satellite->flight_plan().number_of_manœuvres());
  EXPECT_EQ(1 * Newton, satellite->flight_plan().GetManœuvre(0).thrust());
  plugin_->ForgetAllHistoriesBefore(HistoryTime(time, 5));
  EXPECT_LE(HistoryTime(time, 5), satellite->flight_plan().initial_time());
  EXPECT_LE(HistoryTime(time, 5), satellite->history().Begin().time());
  EXPECT_EQ(0, satellite->flight_plan().number_of_manœuvres());
}

TEST_F(PluginTest, ForgetAllHistoriesBeforeAfterPredictionFork) {
  GUID const guid = "Test Satellite";

  InsertAllSolarSystemBodies();
  EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
      .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
  plugin_->EndInitialization();

  EXPECT_CALL(plugin_->mock_ephemeris(), t_max())
      .WillRepeatedly(Return(Instant()));
  EXPECT_CALL(plugin_->mock_ephemeris(), empty()).WillRepeatedly(Return(false));
  EXPECT_CALL(plugin_->mock_ephemeris(), trajectory(_))
      .WillOnce(Return(plugin_->trajectory(SolarSystemFactory::Sun)));
  EXPECT_CALL(plugin_->mock_ephemeris(), Prolong(_)).Times(AnyNumber());
  EXPECT_CALL(plugin_->mock_ephemeris(), FlowWithAdaptiveStep(_, _, _, _, _))
      .WillRepeatedly(DoAll(AppendToDiscreteTrajectory(), Return(true)));
  EXPECT_CALL(plugin_->mock_ephemeris(), FlowWithFixedStep(_, _, _, _))
      .WillRepeatedly(AppendToDiscreteTrajectories());
  EXPECT_CALL(plugin_->mock_ephemeris(), planetary_integrator())
      .WillRepeatedly(
          ReturnRef(McLachlanAtela1992Order5Optimal<Position<Barycentric>>()));

  plugin_->SetPlottingFrame(plugin_->NewBodyCentredNonRotatingNavigationFrame(
      SolarSystemFactory::Sun));
  plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::Earth);
  plugin_->SetVesselStateOffset(guid,
                                RelativeDegreesOfFreedom<AliceSun>(
                                    satellite_initial_displacement_,
                                    satellite_initial_velocity_));

  Instant const& time = initial_time_ + 1 * Second;
  EXPECT_CALL(plugin_->mock_ephemeris(), ForgetBefore(HistoryTime(time, 5)))
      .Times(1);
  plugin_->AdvanceTime(time, Angle());
  plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::Earth);
  plugin_->AdvanceTime(HistoryTime(time, 3), Angle());
  plugin_->UpdatePrediction(guid);
  plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::Earth);
  plugin_->AdvanceTime(HistoryTime(time, 6), Angle());
  plugin_->ForgetAllHistoriesBefore(HistoryTime(time, 5));
  auto const rendered_prediction =
      plugin_->RenderedPrediction(guid, World::origin);
}

TEST_F(PluginDeathTest, VesselFromParentError) {
  GUID const guid = "Test Satellite";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->VesselFromParent(guid);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
        .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
    plugin_->EndInitialization();
    plugin_->VesselFromParent(guid);
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
        .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, SolarSystemFactory::Sun);
    plugin_->VesselFromParent(guid);
  }, "not given an initial state");
}

TEST_F(PluginDeathTest, CelestialFromParentError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->CelestialFromParent(SolarSystemFactory::Earth);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
        .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
    plugin_->EndInitialization();
    plugin_->CelestialFromParent(not_a_body);
  }, "Map key not found");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
        .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
    plugin_->EndInitialization();
    plugin_->CelestialFromParent(SolarSystemFactory::Sun);
  }, "is the sun");
}

TEST_F(PluginTest, VesselInsertionAtInitialization) {
  GUID const guid = "Test Satellite";
  InsertAllSolarSystemBodies();
  EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
      .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
  plugin_->EndInitialization();
  bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                    SolarSystemFactory::Earth);
  EXPECT_TRUE(inserted);
  EXPECT_CALL(plugin_->mock_ephemeris(), Prolong(initial_time_))
      .Times(AnyNumber());
  plugin_->SetVesselStateOffset(guid,
                                RelativeDegreesOfFreedom<AliceSun>(
                                    satellite_initial_displacement_,
                                    satellite_initial_velocity_));
  EXPECT_THAT(plugin_->VesselFromParent(guid),
              Componentwise(
                  AlmostEquals(satellite_initial_displacement_, 13556),
                  AlmostEquals(satellite_initial_velocity_, 1)));
}

TEST_F(PluginTest, UpdateCelestialHierarchy) {
  InsertAllSolarSystemBodies();
  EXPECT_CALL(plugin_->mock_ephemeris(), WriteToMessage(_))
      .WillOnce(SetArgPointee<0>(valid_ephemeris_message_));
  plugin_->EndInitialization();
  EXPECT_CALL(plugin_->mock_ephemeris(), Prolong(_)).Times(AnyNumber());
  for (int index = SolarSystemFactory::Sun + 1;
       index <= SolarSystemFactory::LastMajorBody;
       ++index) {
    plugin_->UpdateCelestialHierarchy(index, SolarSystemFactory::Sun);
  }
  for (int index = SolarSystemFactory::Sun + 1;
       index <= SolarSystemFactory::LastMajorBody;
       ++index) {
    auto const to_icrf = id_icrf_barycentric_.orthogonal_map().Inverse() *
                         plugin_->InversePlanetariumRotation().Forget();
    RelativeDegreesOfFreedom<ICRFJ2000Equator> const from_parent =
        solar_system_->initial_state(SolarSystemFactory::name(index)) -
        solar_system_->initial_state(
            SolarSystemFactory::name(SolarSystemFactory::Sun));
    // All these worlds are fine -- except Triton.
    // Attempt no computation there.
    if (index == SolarSystemFactory::Triton) {
      EXPECT_THAT(
          from_parent,
          Componentwise(
              AlmostEquals(
                  to_icrf(plugin_->CelestialFromParent(index).displacement()),
                  2),
              AlmostEquals(
                  to_icrf(plugin_->CelestialFromParent(index).velocity()),
                  20159592)))
          << SolarSystemFactory::name(index);
    } else {
      // TODO(egg): I'm not sure that's really fine, this seems to be 20 bits?
      // and 24 on Triton.  What's going on?
      EXPECT_THAT(
          from_parent,
          Componentwise(
              AlmostEquals(
                  to_icrf(plugin_->CelestialFromParent(index).displacement()),
                  0, 50),
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
                initial_time_,
                0 * Radian);
  plugin.InsertCelestialJacobiKeplerian(
      SolarSystemFactory::Sun,
      /*parent_index=*/std::experimental::nullopt,
      /*keplerian_elements=*/std::experimental::nullopt,
      std::move(sun_body_));
  plugin.EndInitialization();
  not_null<std::unique_ptr<NavigationFrame>> navigation_frame =
      plugin.NewBodyCentredNonRotatingNavigationFrame(SolarSystemFactory::Sun);
  not_null<const NavigationFrame*> const navigation_frame_copy =
      navigation_frame.get();
  plugin.SetPlottingFrame(std::move(navigation_frame));
  EXPECT_EQ(navigation_frame_copy, plugin.GetPlottingFrame());
  Vector<double, Navball> x_navball({1, 0, 0});
  Vector<double, Navball> y_navball({0, 1, 0});
  Vector<double, Navball> z_navball({0, 0, 1});
  Vector<double, World> x_world({-1, 0, 0});
  Vector<double, World> y_world({0, 1, 0});
  Vector<double, World> z_world({0, 0, -1});
  auto const navball = plugin.NavballFrameField(World::origin);
  EXPECT_THAT(
      AbsoluteError(x_world, navball->FromThisFrame(World::origin)(x_navball)),
      VanishesBefore(1, 4));
  EXPECT_THAT(
      AbsoluteError(y_world, navball->FromThisFrame(World::origin)(y_navball)),
      VanishesBefore(1, 0));
  EXPECT_THAT(
      AbsoluteError(z_world, navball->FromThisFrame(World::origin)(z_navball)),
      VanishesBefore(1, 4));
}

TEST_F(PluginTest, Frenet) {
  // Create a plugin with planetarium rotation 0.
  Plugin plugin(initial_time_,
                initial_time_,
                0 * Radian);
  auto sun_body = make_not_null_unique<RotatingBody<Barycentric>>(
      MassiveBody::Parameters(solar_system_->gravitational_parameter(
          SolarSystemFactory::name(SolarSystemFactory::Earth))),
      body_rotation_);
  plugin.InsertCelestialJacobiKeplerian(
      SolarSystemFactory::Earth,
      /*parent_index=*/std::experimental::nullopt,
      /*keplerian_elements=*/std::experimental::nullopt,
      std::move(sun_body));
  plugin.EndInitialization();
  Permutation<AliceSun, World> const alice_sun_to_world =
      Permutation<AliceSun, World>(Permutation<AliceSun, World>::XZY);
  GUID const satellite = "satellite";
  plugin.InsertOrKeepVessel(satellite, SolarSystemFactory::Earth);
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
          SolarSystemFactory::Earth);
  EXPECT_THAT(plugin.VesselTangent(satellite), AlmostEquals(t, 5, 17));
  EXPECT_THAT(plugin.VesselNormal(satellite), AlmostEquals(n, 4, 11));
  EXPECT_THAT(plugin.VesselBinormal(satellite), AlmostEquals(b, 1, 6));
  EXPECT_THAT(
      plugin.VesselVelocity(satellite),
      AlmostEquals(alice_sun_to_world(satellite_initial_velocity_), 7, 15));
}

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
