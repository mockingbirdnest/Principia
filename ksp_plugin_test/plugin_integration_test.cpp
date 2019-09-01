
#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "astronomy/frames.hpp"
#include "astronomy/time_scales.hpp"
#include "geometry/identity.hpp"
#include "geometry/permutation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "physics/massive_body.hpp"
#include "quantities/astronomy.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

using astronomy::ICRS;
using astronomy::ParseTT;
using base::make_not_null_unique;
using geometry::AffineMap;
using geometry::Bivector;
using geometry::Identity;
using geometry::Permutation;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using physics::KeplerianElements;
using physics::MassiveBody;
using physics::SolarSystem;
using quantities::Abs;
using quantities::Acceleration;
using quantities::ArcTan;
using quantities::Cos;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::NaN;
using quantities::Pow;
using quantities::Sin;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Sqrt;
using quantities::astronomy::AstronomicalUnit;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Minute;
using quantities::si::Radian;
using testing_utilities::AbsoluteError;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystemFactory;
using testing_utilities::operator""_⑴;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::Le;
using ::testing::Lt;

namespace {

PartId const part_id = 789;
GUID const vessel_guid = "123-456";
std::string const part_name = "Picard's desk";
std::string const vessel_name = "NCC-1701-D";

}  // namespace

class PluginIntegrationTest : public testing::Test {
 protected:
  PluginIntegrationTest()
      : looking_glass_(Permutation<ICRS, AliceSun>::XZY),
        solar_system_(
            SolarSystemFactory::AtСпутник1Launch(
                SolarSystemFactory::Accuracy::MinorAndMajorBodies)),
        initial_time_("JD2451545.0625"),
        planetarium_rotation_(1 * Radian),
        plugin_(make_not_null_unique<Plugin>(initial_time_,
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
  }

  void InsertAllSolarSystemBodies() {
    for (int index = SolarSystemFactory::Sun;
         index <= SolarSystemFactory::LastBody;
         ++index) {
      std::optional<Index> parent_index =
          index == SolarSystemFactory::Sun
              ? std::nullopt
              : std::make_optional(SolarSystemFactory::parent(index));
      plugin_->InsertCelestialAbsoluteCartesian(
          index,
          parent_index,
          solar_system_->gravity_model_message(
              SolarSystemFactory::name(index)),
          solar_system_->cartesian_initial_state_message(
              SolarSystemFactory::name(index)));
    }
  }

  Permutation<ICRS, AliceSun> looking_glass_;
  not_null<std::unique_ptr<SolarSystem<ICRS>>> solar_system_;
  std::string initial_time_;
  Angle planetarium_rotation_;

  not_null<std::unique_ptr<Plugin>> plugin_;

  // These initial conditions will yield a low circular orbit around Earth.
  Displacement<AliceSun> satellite_initial_displacement_;
  Velocity<AliceSun> satellite_initial_velocity_;
};

TEST_F(PluginIntegrationTest, AdvanceTimeWithCelestialsOnly) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
#if defined(_DEBUG)
  Time const δt = 2 * Second;
#else
  Time const δt = 0.02 * Second;
#endif
  Angle const planetarium_rotation = 42 * Radian;
  // We step for long enough that we will find a new segment.
  Instant const initial_time = ParseTT(initial_time_);
  Instant t = initial_time;
  for (t += δt; t < initial_time + 10 * 45 * Minute; t += δt) {
    plugin_->AdvanceTime(t, planetarium_rotation);
  }
  EXPECT_THAT(
      RelativeError(
          plugin_->CelestialFromParent(
              SolarSystemFactory::Earth).displacement().Norm(),
          1 * AstronomicalUnit),
      Lt(0.01));
  serialization::Plugin plugin_message;
  plugin_->WriteToMessage(&plugin_message);
  plugin_ = Plugin::ReadFromMessage(plugin_message);
  // Having saved and loaded, we compute a new segment again, this probably
  // exercises apocalypse-type bugs.
  for (; t < initial_time + 20 * 45 * Minute; t += δt) {
    plugin_->AdvanceTime(t, planetarium_rotation);
  }
  EXPECT_THAT(
      RelativeError(
          plugin_->CelestialFromParent(
              SolarSystemFactory::Earth).displacement().Norm(),
          1 * AstronomicalUnit),
      Lt(0.01));
}

TEST_F(PluginIntegrationTest, BodyCentredNonrotatingNavigationIntegration) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();

  bool inserted;
  plugin_->InsertOrKeepVessel(vessel_guid,
                              vessel_name,
                              SolarSystemFactory::Earth,
                              /*loaded=*/false,
                              inserted);
  plugin_->InsertUnloadedPart(
      part_id,
      part_name,
      vessel_guid,
      RelativeDegreesOfFreedom<AliceSun>(satellite_initial_displacement_,
                                         satellite_initial_velocity_));
  plugin_->PrepareToReportCollisions();
  plugin_->FreeVesselsAndPartsAndCollectPileUps(20 * Milli(Second));

  plugin_->renderer().SetPlottingFrame(
      plugin_->NewBodyCentredNonRotatingNavigationFrame(
          SolarSystemFactory::Earth));
  // We'll check that our orbit is rendered as circular (actually, we only check
  // that it is rendered within a thin spherical shell around the Earth).
  Length perigee = std::numeric_limits<double>::infinity() * Metre;
  Length apogee = -std::numeric_limits<double>::infinity() * Metre;
  Permutation<AliceSun, World> const alice_sun_to_world =
      Permutation<AliceSun, World>(Permutation<AliceSun, World>::XZY);
  Time const δt_long = 10 * Minute;
#if defined(_DEBUG)
  Time const δt_short = 1 * Minute;
#else
  Time const δt_short = 0.02 * Second;
#endif
  Instant const initial_time = ParseTT(initial_time_);
  Instant t = initial_time + δt_short;
  // Exercise #267 by having small time steps at the beginning of the trajectory
  // that are not synchronized with those of the Earth.
  for (; t < initial_time + δt_long; t += δt_short) {
    plugin_->AdvanceTime(
        t,
        1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time));
    plugin_->InsertOrKeepVessel(vessel_guid,
                                vessel_name,
                                SolarSystemFactory::Earth,
                                /*loaded=*/false,
                                inserted);
    VesselSet collided_vessels;
    plugin_->CatchUpLaggingVessels(collided_vessels);
  }
  for (; t < initial_time + 12 * Hour; t += δt_long) {
    plugin_->AdvanceTime(
        t,
        1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time));
    VesselSet collided_vessels;
    plugin_->CatchUpLaggingVessels(collided_vessels);
    plugin_->InsertOrKeepVessel(vessel_guid,
                                vessel_name,
                                SolarSystemFactory::Earth,
                                /*loaded=*/false,
                                inserted);
    // We give the sun an arbitrary nonzero velocity in |World|.
    Position<World> const sun_world_position =
        World::origin + Velocity<World>(
            { 0.1 * AstronomicalUnit / Hour,
             -1.0 * AstronomicalUnit / Hour,
              0.0 * AstronomicalUnit / Hour}) * (t - initial_time);
    auto const& psychohistory =
        plugin_->GetVessel(vessel_guid)->psychohistory();
    auto const rendered_trajectory =
        plugin_->renderer().RenderBarycentricTrajectoryInWorld(
            plugin_->CurrentTime(),
            psychohistory.begin(),
            psychohistory.end(),
            sun_world_position,
            plugin_->PlanetariumRotation());
    Position<World> const earth_world_position =
        sun_world_position + alice_sun_to_world(plugin_->CelestialFromParent(
                                 SolarSystemFactory::Earth).displacement());
    for (auto const& [time, degrees_of_freedom] : *rendered_trajectory) {
      Length const distance =
          (degrees_of_freedom.position() - earth_world_position).Norm();
      perigee = std::min(perigee, distance);
      apogee = std::max(apogee, distance);
    }
    EXPECT_THAT(Abs(apogee - perigee), Lt(3 * Metre));
  }
}

TEST_F(PluginIntegrationTest, BarycentricRotatingNavigationIntegration) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();

  bool inserted;
  plugin_->InsertOrKeepVessel(vessel_guid,
                              vessel_name,
                              SolarSystemFactory::Earth,
                              /*loaded=*/false,
                              inserted);
  // A vessel at the Lagrange point L₅.
  RelativeDegreesOfFreedom<AliceSun> const from_the_earth_to_the_moon =
      plugin_->CelestialFromParent(SolarSystemFactory::Moon);
  Displacement<AliceSun> const from_the_earth_to_l5 =
      from_the_earth_to_the_moon.displacement() / 2 -
          Normalize(from_the_earth_to_the_moon.velocity()) *
              from_the_earth_to_the_moon.displacement().Norm() * Sqrt(3) / 2;
  Velocity<AliceSun> const initial_velocity =
      Rotation<AliceSun, AliceSun>(
          π / 3 * Radian,
          Wedge(from_the_earth_to_the_moon.velocity(),
                from_the_earth_to_the_moon.displacement()))(
              from_the_earth_to_the_moon.velocity());
  plugin_->InsertUnloadedPart(part_id,
                              part_name,
                              vessel_guid,
                              {from_the_earth_to_l5, initial_velocity});
  plugin_->PrepareToReportCollisions();
  plugin_->FreeVesselsAndPartsAndCollectPileUps(20 * Milli(Second));

  plugin_->renderer().SetPlottingFrame(
      plugin_->NewBarycentricRotatingNavigationFrame(SolarSystemFactory::Earth,
                                                     SolarSystemFactory::Moon));
  Permutation<AliceSun, World> const alice_sun_to_world =
      Permutation<AliceSun, World>(Permutation<AliceSun, World>::XZY);
  Time const δt_long = 1 * Hour;
#if defined(_DEBUG)
  Time const duration = 12 * Hour;
  Time const δt_short = 20 * Second;
#else
  Time const duration = 20 * Day;
  Time const δt_short = 0.02 * Second;
#endif
  Instant const initial_time = ParseTT(initial_time_);
  Instant t = initial_time + δt_short;
  // Exercise #267 by having small time steps at the beginning of the trajectory
  // that are not synchronized with those of the Earth and Moon.
  for (; t < initial_time + δt_long; t += δt_short) {
    plugin_->AdvanceTime(
        t,
        1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time));
    VesselSet collided_vessels;
    plugin_->CatchUpLaggingVessels(collided_vessels);
    plugin_->InsertOrKeepVessel(vessel_guid,
                                vessel_name,
                                SolarSystemFactory::Earth,
                                /*loaded=*/false,
                                inserted);
  }
  for (; t < initial_time + duration; t += δt_long) {
    plugin_->AdvanceTime(
        t,
        1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time));
    VesselSet collided_vessels;
    plugin_->CatchUpLaggingVessels(collided_vessels);
    plugin_->InsertOrKeepVessel(vessel_guid,
                                vessel_name,
                                SolarSystemFactory::Earth,
                                /*loaded=*/false,
                                inserted);
  }
  plugin_->AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time));
  VesselSet collided_vessels;
  plugin_->CatchUpLaggingVessels(collided_vessels);
  plugin_->InsertOrKeepVessel(vessel_guid,
                              vessel_name,
                              SolarSystemFactory::Earth,
                              /*loaded=*/false,
                              inserted);
  // We give the sun an arbitrary nonzero velocity in |World|.
  Position<World> const sun_world_position =
      World::origin + Velocity<World>(
          { 0.1 * AstronomicalUnit / Hour,
           -1.0 * AstronomicalUnit / Hour,
            0.0 * AstronomicalUnit / Hour}) * (t - initial_time);
    auto const& psychohistory =
        plugin_->GetVessel(vessel_guid)->psychohistory();
    auto const rendered_trajectory =
        plugin_->renderer().RenderBarycentricTrajectoryInWorld(
            plugin_->CurrentTime(),
            psychohistory.begin(),
            psychohistory.end(),
            sun_world_position,
            plugin_->PlanetariumRotation());
    Position<World> const earth_world_position =
        sun_world_position +
        alice_sun_to_world(
            plugin_->CelestialFromParent(SolarSystemFactory::Earth)
                .displacement());
    Position<World> const moon_world_position =
        earth_world_position +
        alice_sun_to_world(
            plugin_->CelestialFromParent(SolarSystemFactory::Moon)
                .displacement());
    Length const earth_moon =
        (moon_world_position - earth_world_position).Norm();
    for (auto const& [time, degrees_of_freedom] : *rendered_trajectory) {
      Position<World> const position = degrees_of_freedom.position();
      Length const satellite_earth = (position - earth_world_position).Norm();
      Length const satellite_moon = (position - moon_world_position).Norm();
      EXPECT_THAT(RelativeError(earth_moon, satellite_earth), Lt(0.0907));
      EXPECT_THAT(RelativeError(earth_moon, satellite_moon), Lt(0.131));
      EXPECT_THAT(RelativeError(satellite_moon, satellite_earth), Lt(0.148));
  }
  // Check that there are no spikes in the rendered trajectory, i.e., that three
  // consecutive points form a sufficiently flat triangle.  This tests issue
  // #256.
  auto it0 = rendered_trajectory->begin();
  CHECK(it0 != rendered_trajectory->end());
  auto it1 = it0;
  ++it1;
  CHECK(it1 != rendered_trajectory->end());
  auto it2 = it1;
  ++it2;
  while (it2 != rendered_trajectory->end()) {
    EXPECT_THAT((it0.degrees_of_freedom().position() -
                 it2.degrees_of_freedom().position())
                    .Norm(),
                Gt(((it0.degrees_of_freedom().position() -
                     it1.degrees_of_freedom().position())
                        .Norm() +
                    (it1.degrees_of_freedom().position() -
                     it2.degrees_of_freedom().position())
                        .Norm()) /
                   1.5))
        << it0.time();
    ++it0;
    ++it1;
    ++it2;
  }
}

#if 0
// The Enterprise D is a low orbit around a massive body with unit gravitational
// parameter, enters the physics bubble, separates, the saucer section reverses
// the direction of its orbit, the physics bubble ends, the two sections meet
// again on the other side of the body, the main section matches its velocity
// with that of the saucer, they are reunited, the physics bubble ends again.
TEST_F(PluginIntegrationTest, PhysicsBubble) {
  GUID const enterprise_d = "NCC-1701-D";
  GUID const enterprise_d_saucer = "NCC-1701-D (saucer)";
  PartId const engineering_section = 0;
  PartId const saucer_section = 1;
  Index const celestial = 0;
  // We use km-day as our unit system because we need the orbit duration to
  // be much larger than 10 s, the fixed step of the histories.
  Time const period = 2 * π * Day;
  double const ε = 1e-10;
  Time const δt = period * ε;
  Length const a = 1 * Kilo(Metre);
  Speed const v0 = 1 * Kilo(Metre) / Day;
  Instant t;
  Plugin plugin(t, t, 0 * Radian);
  auto sun_body = make_not_null_unique<RotatingBody<Barycentric>>(
      MassiveBody::Parameters(1 * Pow<3>(Kilo(Metre)) / Pow<2>(Day)),
      RotatingBody<Barycentric>::Parameters(
          /*mean_radius=*/1 * Metre,
          /*reference_angle=*/1 * Radian,
          /*reference_instant=*/astronomy::J2000,
          /*angular_frequency=*/1 * Radian / Second,
          /*right_ascension_of_pole=*/0 * Degree,
          /*declination_of_pole=*/90 * Degree));
  plugin.InsertCelestialJacobiKeplerian(
      celestial,
      /*parent_index=*/std::nullopt,
      /*keplerian_elements=*/std::nullopt,
      std::move(sun_body));
  plugin.EndInitialization();

  // Step 1: insert the Enterprise.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.InsertUnloadedPart(
      enterprise_d,
      {Displacement<AliceSun>({a, 0 * a, 0 * a}),
       Velocity<AliceSun>({0 * v0, v0, 0 * v0})});
  plugin.AdvanceTime(t, 0 * Radian);

  // Step 2: physics bubble starts.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        engineering_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin + Displacement<World>({a, 0 * a, 0 * a}),
                Velocity<World>({0 * v0, 0 * v0, v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin + Displacement<World>({a, 0 * a, 0 * a}),
                Velocity<World>({0 * v0, 0 * v0, v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
  }
  plugin.AdvanceTime(t, 0 * Radian);

  // Step 3: separation and saucer burn.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.InsertOrKeepVessel(enterprise_d_saucer, celestial);
  // The value of the offset here should be irrelevant, make sure we notice it
  // if it has an influence.
  plugin.InsertUnloadedPart(
      enterprise_d_saucer,
      {Displacement<AliceSun>({NaN<Length>(),
                               NaN<Length>(),
                               NaN<Length>()}),
       Velocity<AliceSun>({NaN<Speed>(),
                           NaN<Speed>(),
                           NaN<Speed>()})});
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        engineering_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin + Displacement<World>({a, 0 * a, 0 * a}),
                Velocity<World>({0 * v0, 0 * v0, v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
  }
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin + Displacement<World>({a, 0 * a, 0 * a}),
                Velocity<World>({0 * v0, 0 * v0, -v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
  }
  plugin.AdvanceTime(t, 0 * Radian);

  // Step 4: end of physics bubble.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.InsertOrKeepVessel(enterprise_d_saucer, celestial);
  plugin.AdvanceTime(t, 0 * Radian);
  EXPECT_THAT(
      RelativeError(Displacement<AliceSun>({a, 0 * a, 0 * a}),
                    plugin.VesselFromParent(enterprise_d).displacement()),
      Lt(100 * ε));
  EXPECT_THAT(RelativeError(
                  Displacement<AliceSun>({a, 0 * a, 0 * a}),
                  plugin.VesselFromParent(enterprise_d_saucer).displacement()),
              Lt(100 * ε));
  EXPECT_THAT(RelativeError(Velocity<AliceSun>({0 * v0, v0, 0 * v0}),
                            plugin.VesselFromParent(enterprise_d).velocity()),
              Lt(100 * ε));
  EXPECT_THAT(
      RelativeError(Velocity<AliceSun>({0 * v0, -v0, 0 * v0}),
                    plugin.VesselFromParent(enterprise_d_saucer).velocity()),
      Lt(100 * ε));

  // Step 5: coming together on the other side.
  t += 0.5 * period;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.InsertOrKeepVessel(enterprise_d_saucer, celestial);
  plugin.AdvanceTime(t, 0 * Radian);
  EXPECT_THAT(
      RelativeError(Displacement<AliceSun>({-a, 0 * a, 0 * a}),
                    plugin.VesselFromParent(enterprise_d).displacement()),
      Lt(100 * ε));
  EXPECT_THAT(RelativeError(
                  Displacement<AliceSun>({-a, 0 * a, 0 * a}),
                  plugin.VesselFromParent(enterprise_d_saucer).displacement()),
              Lt(100 * ε));
  EXPECT_THAT(RelativeError(Velocity<AliceSun>({0 * v0, -v0, 0 * v0}),
                            plugin.VesselFromParent(enterprise_d).velocity()),
              Lt(100 * ε));
  EXPECT_THAT(
      RelativeError(Velocity<AliceSun>({0 * v0, v0, 0 * v0}),
                    plugin.VesselFromParent(enterprise_d_saucer).velocity()),
      Lt(100 * ε));

  // Step 6: reopen physics bubble.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.InsertOrKeepVessel(enterprise_d_saucer, celestial);
  // The absolute world positions don't matter, at least one vessel (indeed all)
  // are pre-existing.  Exercise this.
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        engineering_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin + Displacement<World>({1729 * a, 0 * a, 0 * a}),
                Velocity<World>({0 * v0, 0 * v0, -v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
  }
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin + Displacement<World>({1729 * a, 0 * a, 0 * a}),
                Velocity<World>({0 * v0, 0 * v0, v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
  }
  plugin.AdvanceTime(t, 0 * Radian);

  // Step 7: match velocities.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.InsertOrKeepVessel(enterprise_d_saucer, celestial);
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        engineering_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin, Velocity<World>({0 * v0, 0 * v0, v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
  }
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin, Velocity<World>({0 * v0, 0 * v0, v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
  }
  plugin.AdvanceTime(t, 0 * Radian);

  // Step 8: docking.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        engineering_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin, Velocity<World>({0 * v0, 0 * v0, v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                World::origin, Velocity<World>({0 * v0, 0 * v0, v0})),
            1 * Kilogram, Vector<Acceleration, World>()));
  }
  plugin.AdvanceTime(t, 0 * Radian);

  // Step 9: close physics bubble.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.AdvanceTime(t, 0 * Radian);
  EXPECT_THAT(
      RelativeError(Displacement<AliceSun>({-a, 0 * a, 0 * a}),
                    plugin.VesselFromParent(enterprise_d).displacement()),
      Lt(100 * ε));
  EXPECT_THAT(RelativeError(Velocity<AliceSun>({0 * v0, v0, 0 * v0}),
                            plugin.VesselFromParent(enterprise_d).velocity()),
              Lt(100 * ε));

  // Step 10: orbit a bit.
  t += period;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.AdvanceTime(t, 0 * Radian);
  EXPECT_THAT(
      RelativeError(Displacement<AliceSun>({-a, 0 * a, 0 * a}),
                    plugin.VesselFromParent(enterprise_d).displacement()),
      Lt(100 * ε));
  EXPECT_THAT(RelativeError(Velocity<AliceSun>({0 * v0, v0, 0 * v0}),
                            plugin.VesselFromParent(enterprise_d).velocity()),
              Lt(100 * ε));
}
#endif

// Checks that we correctly predict a full circular orbit around a massive body
// with unit gravitational parameter at unit distance.  Since predictions are
// only computed on |AdvanceTime()|, we advance time by a small amount.
TEST_F(PluginIntegrationTest, Prediction) {
  Index const celestial = 0;
  Plugin plugin("JD2451545.0", "JD2451545.0", 0 * Radian);
  serialization::GravityModel::Body gravity_model;
  CHECK(google::protobuf::TextFormat::ParseFromString(
      R"(name                    : "Celestial"
         gravitational_parameter : "1 m^3/s^2"
         reference_instant       : "JD2451545.0"
         mean_radius             : "1 m"
         axis_right_ascension    : "0 deg"
         axis_declination        : "90 deg"
         reference_angle         : "1 rad"
         angular_frequency       : "1 rad/s")",
      &gravity_model));
  serialization::InitialState::Keplerian::Body initial_state;
  CHECK(google::protobuf::TextFormat::ParseFromString(
      R"(name : "Celestial")",
      &initial_state));
  plugin.InsertCelestialJacobiKeplerian(
      celestial,
      /*parent_index=*/std::nullopt,
      gravity_model,
      initial_state);
  plugin.EndInitialization();

  bool inserted;
  plugin.InsertOrKeepVessel(vessel_guid,
                            vessel_name,
                            celestial,
                            /*loaded=*/false,
                            inserted);
  plugin.InsertUnloadedPart(
      part_id,
      part_name,
      vessel_guid,
      {Displacement<AliceSun>({1 * Metre, 0 * Metre, 0 * Metre}),
       Velocity<AliceSun>(
           {0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second})});
  plugin.PrepareToReportCollisions();
  plugin.FreeVesselsAndPartsAndCollectPileUps(20 * Milli(Second));

  plugin.renderer().SetPlottingFrame(
      plugin.NewBodyCentredNonRotatingNavigationFrame(celestial));
  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters(
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          DormandالمكاوىPrince1986RKN434FM,
          Position<Barycentric>>(),
      /*max_steps=*/14,
      /*length_integration_tolerance=*/1 * Milli(Metre),
      /*speed_integration_tolerance=*/1 * Milli(Metre) / Second);
  plugin.SetPredictionAdaptiveStepParameters(vessel_guid,
                                             adaptive_step_parameters);
  plugin.AdvanceTime(Instant() + 1e-10 * Second, 0 * Radian);

  // Polling for the integration to happen.
  do {
    plugin.UpdatePrediction(vessel_guid);
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(100ms);
  } while (plugin.GetVessel(vessel_guid)->prediction().Size() != 15);

  auto const& prediction = plugin.GetVessel(vessel_guid)->prediction();
  auto const rendered_prediction =
      plugin.renderer().RenderBarycentricTrajectoryInWorld(
          plugin.CurrentTime(),
          prediction.Fork(),
          prediction.end(),
          World::origin,
          plugin.PlanetariumRotation());
  EXPECT_EQ(15, rendered_prediction->Size());
  int index = 0;
  for (auto it = rendered_prediction->begin();
       it != rendered_prediction->end();
       ++it, ++index) {
    auto const& position = it.degrees_of_freedom().position();
    EXPECT_THAT(AbsoluteError((position - World::origin).Norm(), 1 * Metre),
                Lt(0.5 * Milli(Metre)));
    if (index >= 5) {
      EXPECT_THAT(AbsoluteError((position - World::origin).Norm(), 1 * Metre),
                  Gt(0.1 * Milli(Metre)));
    }
  }
  EXPECT_THAT(AbsoluteError(
                  rendered_prediction->rbegin().degrees_of_freedom().position(),
                  Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre}) +
                      World::origin),
              IsNear(29_⑴ * Milli(Metre)));
}

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
