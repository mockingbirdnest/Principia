#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/solar_system.hpp"

namespace principia {

using quantities::Abs;
using quantities::ArcTan;
using quantities::Cos;
using quantities::Sin;
using quantities::Sqrt;
using quantities::si::Day;
using quantities::si::Hour;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::AstronomicalUnit;
using testing_utilities::AbsoluteError;
using testing_utilities::ICRFJ2000Ecliptic;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystem;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::Le;
using ::testing::Lt;

namespace ksp_plugin {

class PluginIntegrationTest : public testing::Test {
 protected:
  PluginIntegrationTest()
      : icrf_to_barycentric_positions_(ICRFJ2000Ecliptic::origin,
                                       Barycentric::origin,
                                       ircf_to_barycentric_linear_),
        looking_glass_(Permutation<ICRFJ2000Ecliptic, AliceSun>::XZY),
        solar_system_(SolarSystem::AtСпутник1Launch(
                      SolarSystem::Accuracy::kAllBodiesAndOblateness)),
        initial_time_(42 * Second),
        planetarium_rotation_(1 * Radian),
        plugin_(make_not_null_unique<Plugin>(initial_time_,
                                             planetarium_rotation_)),
        bodies_(solar_system_->massive_bodies()) {
    sun_gravitational_parameter_ =
        bodies_[SolarSystem::kSun]->gravitational_parameter();
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

  DegreesOfFreedom<Barycentric> ICRFToBarycentric(
      DegreesOfFreedom<ICRFJ2000Ecliptic> const& degrees_of_freedom) {
    return {icrf_to_barycentric_positions_(degrees_of_freedom.position()),
            ircf_to_barycentric_linear_(degrees_of_freedom.velocity())};
  }

  void InsertAllSolarSystemBodies() {
    for (std::size_t index = SolarSystem::kSun;
         index < bodies_.size();
         ++index) {
      std::unique_ptr<Index> parent_index =
          index == SolarSystem::kSun
              ? nullptr
              : std::make_unique<Index>(SolarSystem::parent(index));
      DegreesOfFreedom<Barycentric> const initial_state =
          ICRFToBarycentric(solar_system_->trajectories()[index]->
                            last().degrees_of_freedom());
    // barf_cast, to be fixed when we have rvalue conversion.
    plugin_->DirectlyInsertCelestial(
        index,
        parent_index.get(),
        initial_state,
        std::move(*reinterpret_cast<std::unique_ptr<MassiveBody>*>(
            &bodies_[index])));
    }
  }

  Identity<ICRFJ2000Ecliptic, Barycentric> ircf_to_barycentric_linear_;
  AffineMap<ICRFJ2000Ecliptic,
            Barycentric,
            Length,
            Identity> icrf_to_barycentric_positions_;
  Permutation<ICRFJ2000Ecliptic, AliceSun> looking_glass_;
  not_null<std::unique_ptr<SolarSystem>> solar_system_;
  Instant initial_time_;
  GravitationalParameter sun_gravitational_parameter_;
  Angle planetarium_rotation_;

  not_null<std::unique_ptr<Plugin>> plugin_;

  // These initial conditions will yield a low circular orbit around Earth.
  Displacement<AliceSun> satellite_initial_displacement_;
  Velocity<AliceSun> satellite_initial_velocity_;

 private:
  SolarSystem::Bodies bodies_;
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
  Instant t = initial_time_;
  for (t += δt; t < initial_time_ + 10 * 45 * Minute; t += δt) {
    plugin_->AdvanceTime(t, planetarium_rotation);
  }
  EXPECT_THAT(
      RelativeError(
          plugin_->
              CelestialFromParent(SolarSystem::kEarth).displacement().Norm(),
          1 * AstronomicalUnit),
      Lt(0.01));
  serialization::Plugin plugin_message;
  plugin_->WriteToMessage(&plugin_message);
  plugin_ = Plugin::ReadFromMessage(plugin_message);
  // Having saved and loaded, we compute a new segment again, this probably
  // exercises apocalypse-type bugs.
  for (; t < initial_time_ + 20 * 45 * Minute; t += δt) {
    plugin_->AdvanceTime(t, planetarium_rotation);
  }
  EXPECT_THAT(
      RelativeError(
          plugin_->
              CelestialFromParent(SolarSystem::kEarth).displacement().Norm(),
          1 * AstronomicalUnit),
      Lt(0.01));
}

TEST_F(PluginIntegrationTest, BodyCentredNonrotatingRenderingIntegration) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  GUID const satellite = "satellite";
  plugin_->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  plugin_->SetVesselStateOffset(satellite,
                                RelativeDegreesOfFreedom<AliceSun>(
                                    satellite_initial_displacement_,
                                    satellite_initial_velocity_));
  not_null<std::unique_ptr<RenderingTransforms>> const geocentric =
      plugin_->NewBodyCentredNonRotatingTransforms(SolarSystem::kEarth);
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
  Instant t = initial_time_ + δt_short;
  // Exercise #267 by having small time steps at the beginning of the trajectory
  // that are not synchronized with those of the Earth.
  for (; t < initial_time_ + δt_long; t += δt_short) {
    plugin_->AdvanceTime(
        t,
        1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin_->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  }
  for (; t < initial_time_ + 12 * Hour; t += δt_long) {
    plugin_->AdvanceTime(
        t,
        1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin_->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
    // We give the sun an arbitrary nonzero velocity in |World|.
    Position<World> const sun_world_position =
        World::origin + Velocity<World>(
            { 0.1 * AstronomicalUnit / Hour,
             -1.0 * AstronomicalUnit / Hour,
              0.0 * AstronomicalUnit / Hour}) * (t - initial_time_);
    RenderedTrajectory<World> const rendered_trajectory =
        plugin_->RenderedVesselTrajectory(satellite,
                                          geocentric.get(),
                                          sun_world_position);
    Position<World> const earth_world_position =
        sun_world_position + alice_sun_to_world(
            plugin_->CelestialFromParent(SolarSystem::kEarth).displacement());
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
    EXPECT_THAT(Abs(apogee - perigee), Lt(3 * Metre));
  }
}

TEST_F(PluginIntegrationTest, BarycentricRotatingRenderingIntegration) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  GUID const satellite = "satellite";
  plugin_->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  // A vessel at the Lagrange point L₅.
  RelativeDegreesOfFreedom<AliceSun> const from_the_earth_to_the_moon =
      plugin_->CelestialFromParent(SolarSystem::kMoon);
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
  plugin_->SetVesselStateOffset(satellite,
                                {from_the_earth_to_l5, initial_velocity});
  not_null<std::unique_ptr<RenderingTransforms>> const earth_moon_barycentric =
      plugin_->NewBarycentricRotatingTransforms(SolarSystem::kEarth,
                                                SolarSystem::kMoon);
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
  Instant t = initial_time_ + δt_short;
  // Exercise #267 by having small time steps at the beginning of the trajectory
  // that are not synchronized with those of the Earth and Moon.
  for (; t < initial_time_ + δt_long; t += δt_short) {
    plugin_->AdvanceTime(
        t,
        1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin_->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  }
  for (; t < initial_time_ + duration; t += δt_long) {
    plugin_->AdvanceTime(
        t,
        1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin_->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  }
  plugin_->AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
  plugin_->InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  // We give the sun an arbitrary nonzero velocity in |World|.
  Position<World> const sun_world_position =
      World::origin + Velocity<World>(
          { 0.1 * AstronomicalUnit / Hour,
           -1.0 * AstronomicalUnit / Hour,
            0.0 * AstronomicalUnit / Hour}) * (t - initial_time_);
  RenderedTrajectory<World> const rendered_trajectory =
      plugin_->RenderedVesselTrajectory(satellite,
                                        earth_moon_barycentric.get(),
                                        sun_world_position);
  Position<World> const earth_world_position =
      sun_world_position +
      alice_sun_to_world(
          plugin_->CelestialFromParent(SolarSystem::kEarth).displacement());
  Position<World> const moon_world_position =
      earth_world_position +
      alice_sun_to_world(
          plugin_->CelestialFromParent(SolarSystem::kMoon).displacement());
  Length const earth_moon = (moon_world_position - earth_world_position).Norm();
  for (auto const segment : rendered_trajectory) {
    Length const satellite_earth =
        (segment.begin - earth_world_position).Norm();
    Length const satellite_moon = (segment.begin - moon_world_position).Norm();
    EXPECT_THAT(RelativeError(earth_moon, satellite_earth), Lt(0.0907));
    EXPECT_THAT(RelativeError(earth_moon, satellite_moon), Lt(0.131));
    EXPECT_THAT(RelativeError(satellite_moon, satellite_earth), Lt(0.148));
  }
  // Check continuity.
  for (std::size_t i = 0; i + 1 < rendered_trajectory.size(); ++i) {
    EXPECT_THAT(rendered_trajectory[i].end,
                Eq(rendered_trajectory[i + 1].begin));
  }
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
}

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
  double const ε = 1E-10;
  Time const δt = period * ε;
  Length const a = 1 * Kilo(Metre);
  Speed const v0 = 1 * Kilo(Metre) / Day;
  Instant t;
  Plugin plugin(t, 0 * Radian);
  plugin.InsertSun(celestial, 1 * Pow<3>(Kilo(Metre)) / Pow<2>(Day));
  plugin.EndInitialization();

  // Step 1: insert the Enterprise.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.SetVesselStateOffset(
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
                {World::origin + Displacement<World>({a, 0 * a, 0 * a}),
                 Velocity<World>({0 * v0, 0 * v0, v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                {World::origin + Displacement<World>({a, 0 * a, 0 * a}),
                 Velocity<World>({0 * v0, 0 * v0, v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    plugin.AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
  }
  plugin.AdvanceTime(t, 0 * Radian);
  EXPECT_THAT(plugin.BubbleDisplacementCorrection(World::origin).Norm(),
              Lt(100 * ε * a));
  EXPECT_THAT(plugin.BubbleVelocityCorrection(celestial).Norm(),
              Lt(100 * ε * v0));

  // Step 3: separation and saucer burn.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  plugin.InsertOrKeepVessel(enterprise_d_saucer, celestial);
  // The value of the offset here should be irrelevant, make sure we notice it
  // if it has an influence.
  plugin.SetVesselStateOffset(
      enterprise_d_saucer,
      {Displacement<AliceSun>({std::numeric_limits<double>::quiet_NaN() * a,
                               std::numeric_limits<double>::quiet_NaN() * a,
                               std::numeric_limits<double>::quiet_NaN() * a}),
       Velocity<AliceSun>({std::numeric_limits<double>::quiet_NaN() * v0,
                           std::numeric_limits<double>::quiet_NaN() * v0,
                           std::numeric_limits<double>::quiet_NaN() * v0})});
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        engineering_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                {World::origin + Displacement<World>({a, 0 * a, 0 * a}),
                 Velocity<World>({0 * v0, 0 * v0, v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    plugin.AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
  }
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                {World::origin + Displacement<World>({a, 0 * a, 0 * a}),
                 Velocity<World>({0 * v0, 0 * v0, -v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    plugin.AddVesselToNextPhysicsBubble(enterprise_d_saucer, std::move(parts));
  }
  plugin.AdvanceTime(t, 0 * Radian);
  EXPECT_THAT(plugin.BubbleDisplacementCorrection(World::origin).Norm(),
              Lt(100 * ε * a));
  EXPECT_THAT(plugin.BubbleVelocityCorrection(celestial).Norm(),
              Lt(100 * ε * v0));

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
  // are pre-existing. exercise this.
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        engineering_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                {World::origin + Displacement<World>({1729 * a, 0 * a, 0 * a}),
                 Velocity<World>({0 * v0, 0 * v0, -v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    plugin.AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
  }
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                {World::origin + Displacement<World>({1729 * a, 0 * a, 0 * a}),
                 Velocity<World>({0 * v0, 0 * v0, v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    plugin.AddVesselToNextPhysicsBubble(enterprise_d_saucer, std::move(parts));
  }
  plugin.AdvanceTime(t, 0 * Radian);
  EXPECT_THAT(RelativeError(Displacement<World>({-1730 * a, 0 * a, 0 * a}),
                            plugin.BubbleDisplacementCorrection(World::origin)),
              Lt(100 * ε));
  EXPECT_THAT(plugin.BubbleVelocityCorrection(celestial).Norm(),
              Lt(100 * ε * v0));

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
                {World::origin, Velocity<World>({0 * v0, 0 * v0, v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    plugin.AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
  }
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                {World::origin, Velocity<World>({0 * v0, 0 * v0, v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    plugin.AddVesselToNextPhysicsBubble(enterprise_d_saucer, std::move(parts));
  }
  plugin.AdvanceTime(t, 0 * Radian);
  EXPECT_THAT(
      plugin.BubbleDisplacementCorrection(
          World::origin + Displacement<World>({a, 0 * a, 0 * a})).Norm(),
      Lt(100 * ε * a));
  EXPECT_THAT(plugin.BubbleVelocityCorrection(celestial).Norm(),
              Lt(100 * ε * v0));

  // Step 8: docking.
  t += δt;
  plugin.InsertOrKeepVessel(enterprise_d, celestial);
  {
    std::vector<IdAndOwnedPart> parts;
    parts.emplace_back(
        engineering_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                {World::origin, Velocity<World>({0 * v0, 0 * v0, v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    parts.emplace_back(
        saucer_section,
        make_not_null_unique<Part<World>>(
            DegreesOfFreedom<World>(
                {World::origin, Velocity<World>({0 * v0, 0 * v0, v0})}),
            1 * Kilogram, Vector<Acceleration, World>()));
    plugin.AddVesselToNextPhysicsBubble(enterprise_d, std::move(parts));
  }
  plugin.AdvanceTime(t, 0 * Radian);
  EXPECT_THAT(
      plugin.BubbleDisplacementCorrection(
          World::origin + Displacement<World>({a, 0 * a, 0 * a})).Norm(),
      Lt(100 * ε * a));
  EXPECT_THAT(plugin.BubbleVelocityCorrection(celestial).Norm(),
              Lt(100 * ε * v0));

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

// Checks that we correctly predict a full circular orbit around a massive body
// with unit gravitational parameter at unit distance.  Since predictions are
// only computed on |AdvanceTime()|, we advance time by a small amount.
TEST_F(PluginIntegrationTest, Prediction) {
  GUID const satellite = "satellite";
  Index const celestial = 0;
  Plugin plugin(Instant(), 0 * Radian);
  plugin.InsertSun(celestial, SIUnit<GravitationalParameter>());
  plugin.EndInitialization();
  EXPECT_TRUE(plugin.InsertOrKeepVessel(satellite, celestial));
  auto transforms = plugin.NewBodyCentredNonRotatingTransforms(celestial);
  plugin.SetVesselStateOffset(
      satellite,
      {Displacement<AliceSun>({1 * Metre, 0 * Metre, 0 * Metre}),
       Velocity<AliceSun>(
           {0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second})});
  plugin.set_prediction_length(2 * π * Second);
  plugin.set_prediction_length_tolerance(1 * Milli(Metre));
  plugin.set_prediction_speed_tolerance(1 * Milli(Metre) / Second);
  plugin.AdvanceTime(Instant(1e-10 * Second), 0 * Radian);
  plugin.UpdatePrediction(satellite);
  RenderedTrajectory<World> rendered_prediction =
      plugin.RenderedPrediction(satellite, transforms.get(), World::origin);
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
}

}  // namespace ksp_plugin
}  // namespace principia
