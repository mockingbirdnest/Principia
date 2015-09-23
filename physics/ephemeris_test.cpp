#include "physics/ephemeris.hpp"

#include <limits>
#include <map>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/macros.hpp"
#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/constants.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system.hpp"

namespace principia {

using astronomy::ICRFJ2000Ecliptic;
using astronomy::kSolarSystemBarycentre;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::Abs;
using quantities::ArcTan;
using quantities::Area;
using quantities::Pow;
using quantities::Sqrt;
using quantities::astronomy::LunarDistance;
using quantities::astronomy::SolarMass;
using quantities::constants::GravitationalConstant;
using quantities::si::AstronomicalUnit;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystem;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Ref;

namespace physics {

namespace {

Length const kEarthPolarRadius = 6356.8 * Kilo(Metre);
double const kEarthJ2 = 0.00108262545;

}  // namespace

class EphemerisTest : public testing::Test {
 protected:
  using EarthMoonOrbitPlane = Frame<serialization::Frame::TestTag,
                                    serialization::Frame::TEST, true>;
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  void SetUpEarthMoonSystem(
      not_null<std::vector<not_null<std::unique_ptr<MassiveBody const>>>*> const
          bodies,
      not_null<std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>>*> const
          initial_state,
      not_null<Position<EarthMoonOrbitPlane>*> const centre_of_mass,
      not_null<Time*> const period) {
    // Create the Moon before the Earth to exercise a bug caused by the order of
    // pointers differing from the order of bodies (don't ask).
    auto moon = std::make_unique<MassiveBody const>(7.3459E22 * Kilogram);
    auto earth = std::make_unique<MassiveBody const>(5.9721986E24 * Kilogram);

    // The Earth-Moon system, roughly, with a circular orbit with velocities
    // in the centre-of-mass frame.
    Position<EarthMoonOrbitPlane> const q1(
        Vector<Length, EarthMoonOrbitPlane>({0 * Metre, 0 * Metre, 0 * Metre}));
    Position<EarthMoonOrbitPlane> const q2(
        Vector<Length, EarthMoonOrbitPlane>({0 * Metre,
                                             4E8 * Metre,
                                             0 * Metre}));
    Length const semi_major_axis = (q1 - q2).Norm();
    *period = 2 * π * Sqrt(Pow<3>(semi_major_axis) /
                               (earth->gravitational_parameter() +
                                moon->gravitational_parameter()));
    *centre_of_mass =
        geometry::Barycentre<Vector<Length, EarthMoonOrbitPlane>, Mass>(
            {q1, q2}, {earth->mass(), moon->mass()});
    Velocity<EarthMoonOrbitPlane> const v1(
        {-2 * π * (q1 - *centre_of_mass).Norm() / *period,
         0 * SIUnit<Speed>(),
         0 * SIUnit<Speed>()});
    Velocity<EarthMoonOrbitPlane> const v2(
        {2 * π * (q2 - *centre_of_mass).Norm() / *period,
         0 * SIUnit<Speed>(),
         0 * SIUnit<Speed>()});

    bodies->push_back(std::move(earth));
    bodies->push_back(std::move(moon));
    initial_state->emplace_back(q1, v1);
    initial_state->emplace_back(q2, v2);
  }

  Instant t0_;
};

TEST_F(EphemerisTest, ProlongSpecialCases) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          period / 100,
          5 * Milli(Metre));
  EXPECT_THAT(
      ephemeris.planetary_integrator(),
      Ref(McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>()));

  EXPECT_EQ(t0_ - std::numeric_limits<double>::infinity() * Second,
            ephemeris.t_max());
  EXPECT_EQ(t0_ + std::numeric_limits<double>::infinity() * Second,
            ephemeris.t_min());

  EXPECT_TRUE(ephemeris.empty());
  ephemeris.Prolong(t0_ + period);
  EXPECT_FALSE(ephemeris.empty());
  EXPECT_EQ(t0_, ephemeris.t_min());
  EXPECT_LE(t0_ + period, ephemeris.t_max());
  Instant const t_max = ephemeris.t_max();

  ephemeris.Prolong(t0_ + period / 2);
  EXPECT_EQ(t_max, ephemeris.t_max());

  Instant const last_t =
      geometry::Barycentre<Time, double>({t0_ + period, t_max}, {0.5, 0.5});
  ephemeris.Prolong(last_t);
  EXPECT_EQ(t_max, ephemeris.t_max());
}

// The canonical Earth-Moon system, tuned to produce circular orbits.
TEST_F(EphemerisTest, EarthMoon) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          period / 100,
          5 * Milli(Metre));

  ephemeris.Prolong(t0_ + period);

  ContinuousTrajectory<EarthMoonOrbitPlane> const& earth_trajectory =
      *ephemeris.trajectory(earth);
  ContinuousTrajectory<EarthMoonOrbitPlane> const& moon_trajectory =
      *ephemeris.trajectory(moon);

  ContinuousTrajectory<EarthMoonOrbitPlane>::Hint hint;
  std::vector<Displacement<EarthMoonOrbitPlane>> earth_positions;
  for (int i = 0; i <= 100; ++i) {
    earth_positions.push_back(
      earth_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
          centre_of_mass);
  }
  EXPECT_THAT(earth_positions.size(), Eq(101));
  EXPECT_THAT(Abs(earth_positions[25].coordinates().y), Lt(6E-4 * Metre));
  EXPECT_THAT(Abs(earth_positions[50].coordinates().x), Lt(7E-3 * Metre));
  EXPECT_THAT(Abs(earth_positions[75].coordinates().y), Lt(2E-2 * Metre));
  EXPECT_THAT(Abs(earth_positions[100].coordinates().x), Lt(3E-2 * Metre));

  std::vector<Displacement<EarthMoonOrbitPlane>> moon_positions;
  for (int i = 0; i <= 100; ++i) {
    moon_positions.push_back(
      moon_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
          centre_of_mass);
  }
  EXPECT_THAT(moon_positions.size(), Eq(101));
  EXPECT_THAT(Abs(moon_positions[25].coordinates().y), Lt(5E-2 * Metre));
  EXPECT_THAT(Abs(moon_positions[50].coordinates().x), Lt(6E-1 * Metre));
  EXPECT_THAT(Abs(moon_positions[75].coordinates().y), Lt(2 * Metre));
  EXPECT_THAT(Abs(moon_positions[100].coordinates().x), Lt(2 * Metre));
}

// Test the behavior of ForgetAfter and ForgetBefore on the Earth-Moon system.
TEST_F(EphemerisTest, Forget) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<EarthMoonOrbitPlane>
    ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
      period / 10,
      5 * Milli(Metre));

  ephemeris.Prolong(t0_ + 16 * period);

  ContinuousTrajectory<EarthMoonOrbitPlane> const& earth_trajectory =
    *ephemeris.trajectory(earth);
  ContinuousTrajectory<EarthMoonOrbitPlane> const& moon_trajectory =
    *ephemeris.trajectory(moon);

  Instant t_max = ephemeris.t_max();
  EXPECT_EQ(t0_ + 16 * period, t_max);
  EXPECT_EQ(t_max, earth_trajectory.t_max());
  EXPECT_EQ(t_max, moon_trajectory.t_max());

  ephemeris.ForgetAfter(t0_ + 7 * period);
  t_max = ephemeris.t_max();
  EXPECT_LE(t0_ + 7 * period, t_max);
  EXPECT_GE(t0_ + 7 * period + 180 * Day, t_max);
  EXPECT_EQ(t_max, earth_trajectory.t_max());
  EXPECT_EQ(t_max, moon_trajectory.t_max());

  ephemeris.Prolong(t0_ + 16 * period);
  t_max = ephemeris.t_max();
  EXPECT_EQ(t_max, t0_ + 16 * period);
  EXPECT_EQ(t_max, earth_trajectory.t_max());
  EXPECT_EQ(t_max, moon_trajectory.t_max());

  ephemeris.ForgetAfter(t0_ + 18 * period);
  t_max = ephemeris.t_max();
  EXPECT_EQ(t_max, t0_ + 16 * period);
  EXPECT_EQ(t_max, earth_trajectory.t_max());
  EXPECT_EQ(t_max, moon_trajectory.t_max());

  ephemeris.ForgetBefore(t0_ + 3 * period);
  EXPECT_EQ(t0_ + 3 * period, ephemeris.t_min());
  EXPECT_EQ(t0_ + 3 * period, earth_trajectory.t_min());
  EXPECT_EQ(t0_ + 3 * period, moon_trajectory.t_min());
}

// The Moon alone.  It moves in straight line.
TEST_F(EphemerisTest, Moon) {
  Position<EarthMoonOrbitPlane> const reference_position;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  bodies.erase(bodies.begin());
  initial_state.erase(initial_state.begin());

  MassiveBody const* const moon = bodies[0].get();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          period / 100,
          5 * Milli(Metre));

  ephemeris.Prolong(t0_ + period);

  ContinuousTrajectory<EarthMoonOrbitPlane> const& moon_trajectory =
      *ephemeris.trajectory(moon);

  ContinuousTrajectory<EarthMoonOrbitPlane>::Hint hint;
  DegreesOfFreedom<EarthMoonOrbitPlane> const moon_degrees_of_freedom =
      moon_trajectory.EvaluateDegreesOfFreedom(t0_ + period, &hint);
  Length const q = (moon_degrees_of_freedom.position() -
                    reference_position).coordinates().y;
  Speed const v = moon_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<EarthMoonOrbitPlane>> moon_positions;
  for (int i = 0; i <= 100; ++i) {
    moon_positions.push_back(
        moon_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
            reference_position);
  }

  EXPECT_THAT(moon_positions.size(), Eq(101));
  EXPECT_THAT(moon_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v, 12));
  EXPECT_THAT(moon_positions[25].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v, 11));
  EXPECT_THAT(moon_positions[50].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v, 18));
  EXPECT_THAT(moon_positions[75].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v, 13));
  EXPECT_THAT(moon_positions[100].coordinates().y, Eq(q));
}

// The Earth and a massless probe 1 billion meters away, with the same velocity,
// and an acceleration which exactly compensates gravitational attraction.  Both
// bodies move in straight lines.
TEST_F(EphemerisTest, EarthProbe) {
  Position<EarthMoonOrbitPlane> const reference_position;
  Length const kDistance = 1E9 * Metre;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  bodies.erase(bodies.begin() + 1);
  initial_state.erase(initial_state.begin() + 1);

  MassiveBody const* const earth = bodies[0].get();
  Position<EarthMoonOrbitPlane> const earth_position =
      initial_state[0].position();
  Velocity<EarthMoonOrbitPlane> const earth_velocity =
      initial_state[0].velocity();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          period / 100,
          5 * Milli(Metre));

  MasslessBody probe;
  DiscreteTrajectory<EarthMoonOrbitPlane> trajectory;
  trajectory.Append(t0_,
                    DegreesOfFreedom<EarthMoonOrbitPlane>(
                        earth_position + Vector<Length, EarthMoonOrbitPlane>(
                            {0 * Metre, kDistance, 0 * Metre}),
                        earth_velocity));
  trajectory.set_intrinsic_acceleration(
      [earth, kDistance](Instant const& t) {
    return Vector<Acceleration, EarthMoonOrbitPlane>(
        {0 * SIUnit<Acceleration>(),
         earth->gravitational_parameter() / (kDistance * kDistance),
         0 * SIUnit<Acceleration>()});});

  ephemeris.FlowWithAdaptiveStep(&trajectory,
                                 1E-9 * Metre,
                                 2.6E-15 * Metre / Second,
                                 DormandElMikkawyPrince1986RKN434FM<
                                     Position<EarthMoonOrbitPlane>>(),
                                 t0_ + period);

  ContinuousTrajectory<EarthMoonOrbitPlane> const& earth_trajectory =
      *ephemeris.trajectory(earth);

  ContinuousTrajectory<EarthMoonOrbitPlane>::Hint hint;
  DegreesOfFreedom<EarthMoonOrbitPlane> const earth_degrees_of_freedom =
      earth_trajectory.EvaluateDegreesOfFreedom(t0_ + period, &hint);
  Length const q_earth = (earth_degrees_of_freedom.position() -
                          reference_position).coordinates().y;
  Speed const v_earth = earth_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<EarthMoonOrbitPlane>> earth_positions;
  for (int i = 0; i <= 100; ++i) {
    earth_positions.push_back(
        earth_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
            reference_position);
  }

  EXPECT_THAT(earth_positions.size(), Eq(101));
  EXPECT_THAT(earth_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v_earth, 9));
  EXPECT_THAT(earth_positions[25].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v_earth, 9));
  EXPECT_THAT(earth_positions[50].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v_earth, 4));
  EXPECT_THAT(earth_positions[75].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v_earth, 6));
  EXPECT_THAT(earth_positions[100].coordinates().y, Eq(q_earth));

  Length const q_probe = (trajectory.last().degrees_of_freedom().position() -
                         reference_position).coordinates().y;
  Speed const v_probe =
      trajectory.last().degrees_of_freedom().velocity().coordinates().x;
  std::vector<Displacement<EarthMoonOrbitPlane>> probe_positions;
  for (auto const it : trajectory.Positions()) {
    probe_positions.push_back(it.second - reference_position);
  }
  EXPECT_THAT(probe_positions.size(), Eq(11));
  EXPECT_THAT(probe_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe, 1));
  EXPECT_THAT(probe_positions.back().coordinates().y,
              Eq(q_probe));
}

// The Earth and two massless probes, similar to the previous test but flowing
// with a fixed step.
TEST_F(EphemerisTest, EarthTwoProbes) {
  Position<EarthMoonOrbitPlane> const reference_position;
  Length const kDistance1 = 1E9 * Metre;
  Length const kDistance2 = 3E9 * Metre;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  bodies.erase(bodies.begin() + 1);
  initial_state.erase(initial_state.begin() + 1);

  MassiveBody const* const earth = bodies[0].get();
  Position<EarthMoonOrbitPlane> const earth_position =
      initial_state[0].position();
  Velocity<EarthMoonOrbitPlane> const earth_velocity =
      initial_state[0].velocity();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          period / 100,
          5 * Milli(Metre));

  MasslessBody probe1;
  DiscreteTrajectory<EarthMoonOrbitPlane> trajectory1;
  trajectory1.Append(t0_,
                     DegreesOfFreedom<EarthMoonOrbitPlane>(
                         earth_position + Vector<Length, EarthMoonOrbitPlane>(
                             {0 * Metre, kDistance1, 0 * Metre}),
                         earth_velocity));
  trajectory1.set_intrinsic_acceleration(
      [earth, kDistance1](Instant const& t) {
    return Vector<Acceleration, EarthMoonOrbitPlane>(
        {0 * SIUnit<Acceleration>(),
         earth->gravitational_parameter() / (kDistance1 * kDistance1),
         0 * SIUnit<Acceleration>()});});

  MasslessBody probe2;
  DiscreteTrajectory<EarthMoonOrbitPlane> trajectory2;
  trajectory2.Append(t0_,
                     DegreesOfFreedom<EarthMoonOrbitPlane>(
                         earth_position + Vector<Length, EarthMoonOrbitPlane>(
                             {0 * Metre, -kDistance2, 0 * Metre}),
                         earth_velocity));
  trajectory2.set_intrinsic_acceleration(
      [earth, kDistance2](Instant const& t) {
    return Vector<Acceleration, EarthMoonOrbitPlane>(
        {0 * SIUnit<Acceleration>(),
         -earth->gravitational_parameter() / (kDistance2 * kDistance2),
         0 * SIUnit<Acceleration>()});});

  ephemeris.FlowWithFixedStep({&trajectory1, &trajectory2},
                              period / 1000,
                              t0_ + period);

  ContinuousTrajectory<EarthMoonOrbitPlane> const& earth_trajectory =
      *ephemeris.trajectory(earth);

  ContinuousTrajectory<EarthMoonOrbitPlane>::Hint hint;
  DegreesOfFreedom<EarthMoonOrbitPlane> const earth_degrees_of_freedom =
      earth_trajectory.EvaluateDegreesOfFreedom(t0_ + period, &hint);
  Length const q_earth = (earth_degrees_of_freedom.position() -
                          reference_position).coordinates().y;
  Speed const v_earth = earth_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<EarthMoonOrbitPlane>> earth_positions;
  for (int i = 0; i <= 100; ++i) {
    earth_positions.push_back(
        earth_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
            reference_position);
  }

  EXPECT_THAT(earth_positions.size(), Eq(101));
  EXPECT_THAT(earth_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v_earth, 9));
  EXPECT_THAT(earth_positions[25].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v_earth, 9));
  EXPECT_THAT(earth_positions[50].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v_earth, 4));
  EXPECT_THAT(earth_positions[75].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v_earth, 6));
  EXPECT_THAT(earth_positions[100].coordinates().y, Eq(q_earth));

  Length const q_probe1 = (trajectory1.last().degrees_of_freedom().position() -
                     reference_position).coordinates().y;
  Length const q_probe2 = (trajectory2.last().degrees_of_freedom().position() -
                     reference_position).coordinates().y;
  Speed const v_probe1 =
      trajectory1.last().degrees_of_freedom().velocity().coordinates().x;
  Speed const v_probe2 =
      trajectory2.last().degrees_of_freedom().velocity().coordinates().x;
  std::vector<Displacement<EarthMoonOrbitPlane>> probe1_positions;
  std::vector<Displacement<EarthMoonOrbitPlane>> probe2_positions;
  for (auto const it : trajectory1.Positions()) {
    probe1_positions.push_back(it.second - reference_position);
  }
  for (auto const it : trajectory2.Positions()) {
    probe2_positions.push_back(it.second - reference_position);
  }
#if defined(WE_LOVE_228)
  EXPECT_THAT(probe1_positions.size(), Eq(2));
  EXPECT_THAT(probe2_positions.size(), Eq(2));
#else
  EXPECT_THAT(probe1_positions.size(), Eq(1001));
  EXPECT_THAT(probe2_positions.size(), Eq(1001));
#endif
  EXPECT_THAT(probe1_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe1, 2));
  EXPECT_THAT(probe2_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe2, 2));
  EXPECT_THAT(probe1_positions.back().coordinates().y,
              Eq(q_probe1));
  EXPECT_THAT(probe2_positions.back().coordinates().y,
              Eq(q_probe2));
}

TEST_F(EphemerisTest, Sputnik1ToSputnik2) {
  not_null<std::unique_ptr<SolarSystem>> const at_спутник_1_launch =
      SolarSystem::AtСпутник1Launch(
          SolarSystem::Accuracy::kAllBodiesAndOblateness);
  not_null<std::unique_ptr<SolarSystem>> const at_спутник_2_launch =
      SolarSystem::AtСпутник2Launch(
          SolarSystem::Accuracy::kAllBodiesAndOblateness);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies =
      at_спутник_1_launch->massive_bodies();
  std::vector<DegreesOfFreedom<ICRFJ2000Ecliptic>> const initial_state =
      at_спутник_1_launch->initial_state();

  std::vector<not_null<MassiveBody const*>> unowned_bodies;
  for (auto const& body : bodies) {
    unowned_bodies.push_back(body.get());
  }

  Ephemeris<ICRFJ2000Ecliptic>
      ephemeris(
          std::move(bodies),
          initial_state,
          at_спутник_1_launch->time(),
          McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Ecliptic>>(),
          45 * Minute,
          5 * Milli(Metre));

  ephemeris.Prolong(at_спутник_2_launch->time());

  // Upper bounds, tight to the nearest order of magnitude.
  static std::map<SolarSystem::Index, Angle> const expected_angle_error = {{}};
  static std::map<SolarSystem::Index,
                  double> const expected_parent_distance_error = {{}};
  static std::map<SolarSystem::Index,
                  double> const expected_parent_offset_error = {
      {SolarSystem::kAriel, 1E-3},
      {SolarSystem::kDione, 1E-3},
      {SolarSystem::kIo, 1E-3},
      {SolarSystem::kOberon, 1E-3},
      {SolarSystem::kTethys, 1E-3},
      {SolarSystem::kTitania, 1E-3},
      {SolarSystem::kTriton, 1E-4},
      {SolarSystem::kCharon, 1E-4},
      {SolarSystem::kRhea, 1E-4},
      {SolarSystem::kTitan, 1E-4},
      {SolarSystem::kUmbriel, 1E-4},
      {SolarSystem::kEris, 1E-5},  // NOTE(egg): we may want Dysnomia.
      {SolarSystem::kEuropa, 1E-5},
      {SolarSystem::kGanymede, 1E-5},
      {SolarSystem::kIapetus, 1E-5},
      {SolarSystem::kMoon, 1E-5},  // What is this?
      {SolarSystem::kCallisto, 1E-6},
      {SolarSystem::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystem::kPluto, 1E-6},  // NOTE(egg): We are missing Hydra and Nyx.
      {SolarSystem::kVenus, 1E-7},
      {SolarSystem::kEarth, 1E-8},
      {SolarSystem::kNeptune, 1E-8},
      {SolarSystem::kSaturn, 1E-8},
      {SolarSystem::kUranus, 1E-8},
      {SolarSystem::kJupiter, 1E-9},
      {SolarSystem::kMars, 1E-9}};
  static std::map<SolarSystem::Index, double> const expected_position_error = {
      {SolarSystem::kEris, 1E-5},  // NOTE(egg): we may want Dysnomia.
      {SolarSystem::kCharon, 1E-6},
      {SolarSystem::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystem::kPluto, 1E-6},
      {SolarSystem::kTethys, 1E-6},
      {SolarSystem::kAriel, 1E-7},
      {SolarSystem::kDione, 1E-7},
      {SolarSystem::kIo, 1E-7},
      {SolarSystem::kMoon, 1E-7},
      {SolarSystem::kOberon, 1E-7},
      {SolarSystem::kRhea, 1E-7},
      {SolarSystem::kTitan, 1E-7},
      {SolarSystem::kTitania, 1E-7},
      {SolarSystem::kVenus, 1E-7},
      {SolarSystem::kCallisto, 1E-8},
      {SolarSystem::kEarth, 1E-8},
      {SolarSystem::kEuropa, 1E-8},
      {SolarSystem::kGanymede, 1E-8},
      {SolarSystem::kIapetus, 1E-8},
      {SolarSystem::kNeptune, 1E-8},
      {SolarSystem::kSaturn, 1E-8},
      {SolarSystem::kSun, 1E-8},
      {SolarSystem::kTriton, 1E-8},
      {SolarSystem::kUmbriel, 1E-8},
      {SolarSystem::kUranus, 1E-8},
      {SolarSystem::kJupiter, 1E-9},
      {SolarSystem::kMars, 1E-9}};
  static std::map<SolarSystem::Index, double> const expected_velocity_error = {
      {SolarSystem::kAriel, 1E-3},
      {SolarSystem::kCharon, 1E-3},
      {SolarSystem::kDione, 1E-3},
      {SolarSystem::kIo, 1E-3},
      {SolarSystem::kPluto, 1E-3},
      {SolarSystem::kTethys, 1E-3},
      {SolarSystem::kEuropa, 1E-4},
      {SolarSystem::kOberon, 1E-4},
      {SolarSystem::kRhea, 1E-4},
      {SolarSystem::kTitania, 1E-4},
      {SolarSystem::kTriton, 1E-4},
      {SolarSystem::kUmbriel, 1E-4},
      {SolarSystem::kEris, 1E-5},  // NOTE(egg): we may want Dysnomia.
      {SolarSystem::kGanymede, 1E-5},
      {SolarSystem::kTitan, 1E-5},
      {SolarSystem::kUranus, 1E-5},
      {SolarSystem::kCallisto, 1E-6},
      {SolarSystem::kIapetus, 1E-6},
      {SolarSystem::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystem::kMoon, 1E-6},
      {SolarSystem::kSaturn, 1E-6},
      {SolarSystem::kEarth, 1E-7},
      {SolarSystem::kJupiter, 1E-7},
      {SolarSystem::kNeptune, 1E-7},
      {SolarSystem::kSun, 1E-7},
      {SolarSystem::kVenus, 1E-7},
      {SolarSystem::kMars, 1E-8}};

  std::vector<DegreesOfFreedom<ICRFJ2000Ecliptic>> const final_state =
      at_спутник_2_launch->initial_state();
  for (std::size_t i = 0; i < unowned_bodies.size(); ++i) {
    SolarSystem::Index const index = static_cast<SolarSystem::Index>(i);
    ContinuousTrajectory<ICRFJ2000Ecliptic> const& trajectory =
        *ephemeris.trajectory(unowned_bodies[i]);
    double const position_error = RelativeError(
        final_state[i].position() - kSolarSystemBarycentre,
        trajectory.EvaluatePosition(at_спутник_2_launch->time(), nullptr) -
            kSolarSystemBarycentre);
    double const velocity_error = RelativeError(
        final_state[i].velocity(),
        trajectory.EvaluateVelocity(at_спутник_2_launch->time(), nullptr));
    EXPECT_THAT(position_error, Lt(expected_position_error.at(index)))
        << SolarSystem::name(i);
    EXPECT_THAT(position_error, Gt(expected_position_error.at(index) / 10.0))
        << SolarSystem::name(i);
    EXPECT_THAT(velocity_error, Lt(expected_velocity_error.at(index)))
        << SolarSystem::name(i);
    EXPECT_THAT(velocity_error, Gt(expected_velocity_error.at(index) / 10.0))
        << SolarSystem::name(i);
    if (i != SolarSystem::kSun) {
      // Look at the error in the position relative to the parent.
      Vector<Length, ICRFJ2000Ecliptic> expected =
          at_спутник_2_launch->trajectories()[i]->
              last().degrees_of_freedom().position() -
          at_спутник_2_launch->trajectories()[SolarSystem::parent(i)]->
              last().degrees_of_freedom().position();
      Vector<Length, ICRFJ2000Ecliptic> actual =
          trajectory.EvaluatePosition(at_спутник_2_launch->time(), nullptr) -
          ephemeris.trajectory(unowned_bodies[SolarSystem::parent(i)])->
              EvaluatePosition(at_спутник_2_launch->time(), nullptr);
      if (expected_angle_error.find(index) != expected_angle_error.end()) {
        Area const product_of_norms = expected.Norm() * actual.Norm();
        Angle const angle = ArcTan(
            Wedge(expected, actual).Norm() / product_of_norms,
            InnerProduct(expected, actual) / product_of_norms);
        EXPECT_THAT(angle / Degree,
                    Gt(expected_angle_error.at(index) / Degree * 0.9))
            << SolarSystem::name(i);
        EXPECT_THAT(angle / Degree,
                    Lt(expected_angle_error.at(index) / Degree * 1.1))
            << SolarSystem::name(i);
      }
      if (expected_parent_distance_error.find(index) !=
          expected_parent_distance_error.end()) {
        double const parent_distance_error = RelativeError(expected.Norm(),
                                                  actual.Norm());
        EXPECT_THAT(parent_distance_error,
                    Lt(expected_parent_distance_error.at(index)))
            << SolarSystem::name(i);
        EXPECT_THAT(parent_distance_error,
                    Gt(expected_parent_distance_error.at(index) / 10.0))
            << SolarSystem::name(i);
      }
      if (expected_parent_offset_error.find(index) !=
          expected_parent_offset_error.end()) {
        double const parent_offset_error =  RelativeError(expected, actual);
        EXPECT_THAT(parent_offset_error,
                    Lt(expected_parent_offset_error.at(index)))
            << SolarSystem::name(i);
        EXPECT_THAT(parent_offset_error,
                    Gt(expected_parent_offset_error.at(index) / 10.0))
            << SolarSystem::name(i);
      }
    }
  }
}

TEST_F(EphemerisTest, Serialization) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          period / 100,
          5 * Milli(Metre));
  ephemeris.Prolong(t0_ + period);

  serialization::Ephemeris message;
  ephemeris.WriteToMessage(&message);

  auto const ephemeris_read =
      Ephemeris<EarthMoonOrbitPlane>::ReadFromMessage(message);
  MassiveBody const* const earth_read = ephemeris_read->bodies()[0];
  MassiveBody const* const moon_read = ephemeris_read->bodies()[1];


  EXPECT_EQ(ephemeris.t_min(), ephemeris_read->t_min());
  EXPECT_EQ(ephemeris.t_max(), ephemeris_read->t_max());
  for (Instant time = ephemeris.t_min();
       time <= ephemeris.t_max();
       time += (ephemeris.t_max() - ephemeris.t_min()) / 100) {
    EXPECT_EQ(ephemeris.trajectory(earth)->EvaluateDegreesOfFreedom(
                  time, nullptr /*hint*/),
              ephemeris_read->trajectory(earth_read)->EvaluateDegreesOfFreedom(
                  time, nullptr /*hint*/));
    EXPECT_EQ(ephemeris.trajectory(moon)->EvaluateDegreesOfFreedom(
                  time, nullptr /*hint*/),
              ephemeris_read->trajectory(moon_read)->EvaluateDegreesOfFreedom(
                  time, nullptr /*hint*/));
  }

  serialization::Ephemeris other_message;
  ephemeris_read->WriteToMessage(&other_message);

  EXPECT_EQ(other_message.SerializeAsString(), message.SerializeAsString());
}

// The gravitational acceleration on at elephant located at the pole.
TEST_F(EphemerisTest, ComputeGravitationalAccelerationMasslessBody) {
  Position<EarthMoonOrbitPlane> const reference_position;
  Time const kDuration = 1 * Second;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;

  auto* earth = new OblateBody<EarthMoonOrbitPlane>(
                        5.9721986E24 * Kilogram,
                        RotatingBody<EarthMoonOrbitPlane>::Parameters(
                            1 * Radian,
                            t0_,
                            AngularVelocity<EarthMoonOrbitPlane>({
                                0 * Radian / Second,
                                0 * Radian / Second,
                                4 * Radian / Second})),
                        OblateBody<EarthMoonOrbitPlane>::Parameters(
                            kEarthJ2,
                            kEarthPolarRadius));
  Velocity<EarthMoonOrbitPlane> const v({0 * SIUnit<Speed>(),
                                         0 * SIUnit<Speed>(),
                                         0 * SIUnit<Speed>()});
  Position<EarthMoonOrbitPlane> const q(
      Vector<Length, EarthMoonOrbitPlane>({0 * AstronomicalUnit,
                                           0 * AstronomicalUnit,
                                           0 * AstronomicalUnit}));

  bodies.push_back(std::unique_ptr<MassiveBody const>(earth));
  initial_state.emplace_back(q, v);

  Position<EarthMoonOrbitPlane> const earth_position =
      initial_state[0].position();
  Velocity<EarthMoonOrbitPlane> const earth_velocity =
      initial_state[0].velocity();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          kDuration / 100,
          5 * Milli(Metre));

  MasslessBody elephant;
  DiscreteTrajectory<EarthMoonOrbitPlane> trajectory;
  trajectory.Append(t0_,
                    DegreesOfFreedom<EarthMoonOrbitPlane>(
                        earth_position + Vector<Length, EarthMoonOrbitPlane>(
                            {0 * Metre, 0 * Metre, kEarthPolarRadius}),
                        earth_velocity));

  ephemeris.FlowWithAdaptiveStep(&trajectory,
                                 1E-9 * Metre,
                                 2.6E-15 * Metre / Second,
                                 DormandElMikkawyPrince1986RKN434FM<
                                     Position<EarthMoonOrbitPlane>>(),
                                 t0_ + kDuration);

  Speed const v_elephant =
      trajectory.last().degrees_of_freedom().velocity().coordinates().x;
  std::vector<Displacement<EarthMoonOrbitPlane>> elephant_positions;
  for (auto const pair : trajectory.Positions()) {
    elephant_positions.push_back(pair.second - reference_position);
  }
  EXPECT_THAT(elephant_positions.size(), Eq(8));
  EXPECT_THAT(elephant_positions.back().coordinates().x,
              AlmostEquals(kDuration * v_elephant, 0));
  EXPECT_LT(RelativeError(elephant_positions.back().coordinates().z,
                          kEarthPolarRadius), 8E-7);

  std::vector<Vector<Acceleration, EarthMoonOrbitPlane>> elephant_accelerations;
  for (auto const& t : trajectory.Times()) {
    elephant_accelerations.push_back(
        ephemeris.ComputeGravitationalAcceleration(&trajectory, t));
  }
  EXPECT_THAT(elephant_accelerations.size(), Eq(8));
  EXPECT_THAT(elephant_accelerations.back().coordinates().x,
              AlmostEquals(0 * SIUnit<Acceleration>(), 0));
  EXPECT_THAT(elephant_accelerations.back().coordinates().y,
              AlmostEquals(0 * SIUnit<Acceleration>(), 0));
  EXPECT_LT(RelativeError(elephant_accelerations.back().coordinates().z,
                          -9.832 * SIUnit<Acceleration>()), 4.9E-5);
}

TEST_F(EphemerisTest, ComputeGravitationalAccelerationMassiveBody) {
  Time const kDuration = 1 * Second;
  double const kJ2 = 1E6;
  Length const kRadius = 1 * LunarDistance;

  Mass const m0 = 1 * SolarMass;
  Mass const m1 = 2 * SolarMass;
  Mass const m2 = 3 * SolarMass;
  Mass const m3 = 4 * SolarMass;

  auto const b0 = new OblateBody<World>(m0,
                                        RotatingBody<World>::Parameters(
                                            1 * Radian,
                                            t0_,
                                            AngularVelocity<World>({
                                                0 * Radian / Second,
                                                0 * Radian / Second,
                                                4 * Radian / Second})),
                                        OblateBody<World>::Parameters(
                                            kJ2, kRadius));
  auto const b1 = new MassiveBody(m1);
  auto const b2 = new MassiveBody(m2);
  auto const b3 = new MassiveBody(m3);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<World>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b0));
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b1));
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b2));
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b3));

  Velocity<World> const v({0 * SIUnit<Speed>(),
                           0 * SIUnit<Speed>(),
                           0 * SIUnit<Speed>()});
  Position<World> const q0(
      Vector<Length, World>({0 * AstronomicalUnit,
                             0 * AstronomicalUnit,
                             0 * AstronomicalUnit}));
  Position<World> const q1(
      Vector<Length, World>({1 * AstronomicalUnit,
                             0 * AstronomicalUnit,
                             0 * AstronomicalUnit}));
  Position<World> const q2(
      Vector<Length, World>({1 * AstronomicalUnit,
                             0 * AstronomicalUnit,
                             1 * AstronomicalUnit}));
  Position<World> const q3(
      Vector<Length, World>({0 * AstronomicalUnit,
                             0 * AstronomicalUnit,
                             1 * AstronomicalUnit}));
  initial_state.emplace_back(q0, v);
  initial_state.emplace_back(q1, v);
  initial_state.emplace_back(q2, v);
  initial_state.emplace_back(q3, v);

  Ephemeris<World>
      ephemeris(std::move(bodies),
                initial_state,
                t0_,
                McLachlanAtela1992Order5Optimal<Position<World>>(),
                kDuration / 100,
                5 * Milli(Metre));
  ephemeris.Prolong(t0_ + kDuration);

  Vector<Acceleration, World> actual_acceleration0 =
      ephemeris.ComputeGravitationalAcceleration(b0, t0_);
  Vector<Acceleration, World> expected_acceleration0 =
      GravitationalConstant * (m1 * (q1 - q0) / Pow<3>((q1 - q0).Norm()) +
                               m2 * (q2 - q0) / Pow<3>((q2 - q0).Norm()) +
                               m3 * (q3 - q0) / Pow<3>((q3 - q0).Norm())) +
      Vector<Acceleration, World>(
          {(1.5 * m1 - (9 / Sqrt(512)) * m2) * GravitationalConstant *
               Pow<2>(kRadius) * kJ2 / Pow<4>((q0 - q1).Norm()),
           0 * SIUnit<Acceleration>(),
           (-3 * m3 + (3 / Sqrt(512)) * m2) * GravitationalConstant *
               Pow<2>(kRadius) * kJ2 / Pow<4>((q0 - q1).Norm())});
  EXPECT_THAT(actual_acceleration0,
              AlmostEquals(expected_acceleration0, 0, 5));

  Vector<Acceleration, World> actual_acceleration1 =
      ephemeris.ComputeGravitationalAcceleration(b1, t0_);
  Vector<Acceleration, World> expected_acceleration1 =
      GravitationalConstant * (m0 * (q0 - q1) / Pow<3>((q0 - q1).Norm()) +
                               m2 * (q2 - q1) / Pow<3>((q2 - q1).Norm()) +
                               m3 * (q3 - q1) / Pow<3>((q3 - q1).Norm())) +
      Vector<Acceleration, World>(
          {-1.5 * GravitationalConstant * m0 * Pow<2>(kRadius) * kJ2 /
               Pow<4>((q0 - q1).Norm()),
           0 * SIUnit<Acceleration>(),
           0 * SIUnit<Acceleration>()});
  EXPECT_THAT(actual_acceleration1,
              AlmostEquals(expected_acceleration1, 0, 2));

  Vector<Acceleration, World> actual_acceleration2 =
      ephemeris.ComputeGravitationalAcceleration(b2, t0_);
  Vector<Acceleration, World> expected_acceleration2 =
      GravitationalConstant * (m0 * (q0 - q2) / Pow<3>((q0 - q2).Norm()) +
                               m1 * (q1 - q2) / Pow<3>((q1 - q2).Norm()) +
                               m3 * (q3 - q2) / Pow<3>((q3 - q2).Norm())) +
      Vector<Acceleration, World>(
          {(9 / Sqrt(512)) * GravitationalConstant * m0 *
               Pow<2>(kRadius) * kJ2 / Pow<4>((q0 - q1).Norm()),
           0 * SIUnit<Acceleration>(),
           (-3 / Sqrt(512)) * GravitationalConstant * m0 *
               Pow<2>(kRadius) * kJ2 / Pow<4>((q0 - q1).Norm())});
  EXPECT_THAT(actual_acceleration2,
              AlmostEquals(expected_acceleration2, 0, 2));

  Vector<Acceleration, World> actual_acceleration3 =
      ephemeris.ComputeGravitationalAcceleration(b3, t0_);
  Vector<Acceleration, World> expected_acceleration3 =
      GravitationalConstant * (m0 * (q0 - q3) / Pow<3>((q0 - q3).Norm()) +
                               m1 * (q1 - q3) / Pow<3>((q1 - q3).Norm()) +
                               m2 * (q2 - q3) / Pow<3>((q2 - q3).Norm())) +
      Vector<Acceleration, World>(
          {0 * SIUnit<Acceleration>(),
           0 * SIUnit<Acceleration>(),
           3 * GravitationalConstant * m0 * Pow<2>(kRadius) * kJ2 /
               Pow<4>((q0 - q1).Norm())});
  EXPECT_THAT(actual_acceleration3,
              AlmostEquals(expected_acceleration3, 0, 4));
}

}  // namespace physics
}  // namespace principia
