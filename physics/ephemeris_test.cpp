
#include "physics/ephemeris.hpp"

#include <limits>
#include <map>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/macros.hpp"
#include "geometry/barycentre_calculator.hpp"
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
#include "testing_utilities/solar_system_factory.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using astronomy::kSolarSystemBarycentreEquator;
using geometry::Barycentre;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::Abs;
using quantities::ArcTan;
using quantities::Area;
using quantities::Pow;
using quantities::Sqrt;
using quantities::astronomy::JulianYear;
using quantities::astronomy::LunarDistance;
using quantities::astronomy::SolarMass;
using quantities::constants::GravitationalConstant;
using quantities::si::AstronomicalUnit;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystemFactory;
using testing_utilities::VanishesBefore;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Ref;

namespace physics {

namespace {

Length constexpr kEarthPolarRadius = 6356.8 * Kilo(Metre);
double constexpr kEarthJ2 = 0.00108262545;
int constexpr kMaxSteps = 1000;

}  // namespace

class EphemerisTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  EphemerisTest() {
    solar_system_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2433282_500000000.proto.txt");
    t0_ = solar_system_.epoch();
  }

  void SetUpEarthMoonSystem(
      not_null<std::vector<not_null<std::unique_ptr<MassiveBody const>>>*> const
          bodies,
      not_null<std::vector<DegreesOfFreedom<ICRFJ2000Equator>>*> const
          initial_state,
      not_null<Position<ICRFJ2000Equator>*> const centre_of_mass,
      not_null<Time*> const period) {
    // Make the bodies non-oblate so that the system can be computed explicitly.
    serialization::GravityModel::Body earth_gravity_model =
        solar_system_.gravity_model_message("Earth");
    earth_gravity_model.clear_j2();
    earth_gravity_model.clear_reference_radius();
    serialization::GravityModel::Body moon_gravity_model =
        solar_system_.gravity_model_message("Moon");
    moon_gravity_model.clear_j2();
    moon_gravity_model.clear_reference_radius();

    // Create the Moon before the Earth to exercise a bug caused by the order of
    // pointers differing from the order of bodies (don't ask).
    auto moon =
        SolarSystem<ICRFJ2000Equator>::MakeMassiveBody(moon_gravity_model);
    auto earth =
        SolarSystem<ICRFJ2000Equator>::MakeMassiveBody(earth_gravity_model);

    // The Earth-Moon system, roughly, with a circular orbit with velocities
    // in the centre-of-mass frame.
    Position<ICRFJ2000Equator> const q1 = ICRFJ2000Equator::origin +
        Displacement<ICRFJ2000Equator>({0 * Metre, 0 * Metre, 0 * Metre});
    Position<ICRFJ2000Equator> const q2 = ICRFJ2000Equator::origin +
        Displacement<ICRFJ2000Equator>({0 * Metre,
                                        4E8 * Metre,
                                        0 * Metre});
    Length const semi_major_axis = (q1 - q2).Norm();
    *period = 2 * π * Sqrt(Pow<3>(semi_major_axis) /
                               (earth->gravitational_parameter() +
                                moon->gravitational_parameter()));
    *centre_of_mass =
        Barycentre<Position<ICRFJ2000Equator>, Mass>(
            {q1, q2}, {earth->mass(), moon->mass()});
    Velocity<ICRFJ2000Equator> const v1(
        {-2 * π * (q1 - *centre_of_mass).Norm() / *period,
         0 * SIUnit<Speed>(),
         0 * SIUnit<Speed>()});
    Velocity<ICRFJ2000Equator> const v2(
        {2 * π * (q2 - *centre_of_mass).Norm() / *period,
         0 * SIUnit<Speed>(),
         0 * SIUnit<Speed>()});

    bodies->push_back(std::move(earth));
    bodies->push_back(std::move(moon));
    initial_state->emplace_back(q1, v1);
    initial_state->emplace_back(q2, v2);
  }

  SolarSystem<ICRFJ2000Equator> solar_system_;
  Instant t0_;
};

TEST_F(EphemerisTest, ProlongSpecialCases) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> initial_state;
  Position<ICRFJ2000Equator> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  Ephemeris<ICRFJ2000Equator>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              period / 100));
  EXPECT_THAT(
      ephemeris.planetary_integrator(),
      Ref(McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>()));

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
      Barycentre<Instant, double>({t0_ + period, t_max}, {0.5, 0.5});
  ephemeris.Prolong(last_t);
  EXPECT_EQ(t_max, ephemeris.t_max());
}

// The canonical Earth-Moon system, tuned to produce circular orbits.
TEST_F(EphemerisTest, EarthMoon) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> initial_state;
  Position<ICRFJ2000Equator> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<ICRFJ2000Equator>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              period / 100));

  ephemeris.Prolong(t0_ + period);

  ContinuousTrajectory<ICRFJ2000Equator> const& earth_trajectory =
      *ephemeris.trajectory(earth);
  ContinuousTrajectory<ICRFJ2000Equator> const& moon_trajectory =
      *ephemeris.trajectory(moon);

  ContinuousTrajectory<ICRFJ2000Equator>::Hint hint;
  std::vector<Displacement<ICRFJ2000Equator>> earth_positions;
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

  std::vector<Displacement<ICRFJ2000Equator>> moon_positions;
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
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> initial_state;
  Position<ICRFJ2000Equator> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<ICRFJ2000Equator>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              period / 10));

  ephemeris.Prolong(t0_ + 16 * period);

  ContinuousTrajectory<ICRFJ2000Equator> const& earth_trajectory =
    *ephemeris.trajectory(earth);
  ContinuousTrajectory<ICRFJ2000Equator> const& moon_trajectory =
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
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> initial_state;
  Position<ICRFJ2000Equator> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  bodies.erase(bodies.begin());
  initial_state.erase(initial_state.begin());

  MassiveBody const* const moon = bodies[0].get();

  Ephemeris<ICRFJ2000Equator>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              period / 100));

  ephemeris.Prolong(t0_ + period);

  ContinuousTrajectory<ICRFJ2000Equator> const& moon_trajectory =
      *ephemeris.trajectory(moon);

  ContinuousTrajectory<ICRFJ2000Equator>::Hint hint;
  DegreesOfFreedom<ICRFJ2000Equator> const moon_degrees_of_freedom =
      moon_trajectory.EvaluateDegreesOfFreedom(t0_ + period, &hint);
  Length const q = (moon_degrees_of_freedom.position() -
                    ICRFJ2000Equator::origin).coordinates().y;
  Speed const v = moon_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<ICRFJ2000Equator>> moon_positions;
  for (int i = 0; i <= 100; ++i) {
    moon_positions.push_back(
        moon_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
            ICRFJ2000Equator::origin);
  }

  EXPECT_THAT(moon_positions.size(), Eq(101));
  EXPECT_THAT(moon_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v, 362));
  EXPECT_THAT(moon_positions[25].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v, 135));
  EXPECT_THAT(moon_positions[50].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v, 543));
  EXPECT_THAT(moon_positions[75].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v, 383));
  EXPECT_THAT(moon_positions[100].coordinates().y, Eq(q));
}

// The Earth and a massless probe 1 billion meters away, with the same velocity,
// and an acceleration which exactly compensates gravitational attraction.  Both
// bodies move in straight lines.
TEST_F(EphemerisTest, EarthProbe) {
  Length const kDistance = 1E9 * Metre;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> initial_state;
  Position<ICRFJ2000Equator> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  bodies.erase(bodies.begin() + 1);
  initial_state.erase(initial_state.begin() + 1);

  MassiveBody const* const earth = bodies[0].get();
  Position<ICRFJ2000Equator> const earth_position =
      initial_state[0].position();
  Velocity<ICRFJ2000Equator> const earth_velocity =
      initial_state[0].velocity();

  Ephemeris<ICRFJ2000Equator>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              period / 100));

  MasslessBody probe;
  DiscreteTrajectory<ICRFJ2000Equator> trajectory;
  trajectory.Append(t0_,
                    DegreesOfFreedom<ICRFJ2000Equator>(
                        earth_position + Vector<Length, ICRFJ2000Equator>(
                            {0 * Metre, kDistance, 0 * Metre}),
                        earth_velocity));
  auto const intrinsic_acceleration =
      [earth, kDistance](Instant const& t) {
        return Vector<Acceleration, ICRFJ2000Equator>(
            {0 * SIUnit<Acceleration>(),
             earth->gravitational_parameter() / (kDistance * kDistance),
             0 * SIUnit<Acceleration>()});
      };

  ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      intrinsic_acceleration,
      t0_ + period,
      Ephemeris<ICRFJ2000Equator>::AdaptiveStepParameters(
          DormandElMikkawyPrince1986RKN434FM<Position<ICRFJ2000Equator>>(),
          kMaxSteps,
          1E-9 * Metre,
          2.6E-15 * Metre / Second),
          Ephemeris<ICRFJ2000Equator>::unlimited_max_ephemeris_steps);

  ContinuousTrajectory<ICRFJ2000Equator> const& earth_trajectory =
      *ephemeris.trajectory(earth);

  ContinuousTrajectory<ICRFJ2000Equator>::Hint hint;
  DegreesOfFreedom<ICRFJ2000Equator> const earth_degrees_of_freedom =
      earth_trajectory.EvaluateDegreesOfFreedom(t0_ + period, &hint);
  Length const q_earth = (earth_degrees_of_freedom.position() -
                          ICRFJ2000Equator::origin).coordinates().y;
  Speed const v_earth = earth_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<ICRFJ2000Equator>> earth_positions;
  for (int i = 0; i <= 100; ++i) {
    earth_positions.push_back(
        earth_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
            ICRFJ2000Equator::origin);
  }

  EXPECT_THAT(earth_positions.size(), Eq(101));
  EXPECT_THAT(earth_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v_earth, 551));
  EXPECT_THAT(earth_positions[25].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v_earth, 230));
  EXPECT_THAT(earth_positions[50].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v_earth, 413));
  EXPECT_THAT(earth_positions[75].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v_earth, 622));
  EXPECT_THAT(earth_positions[100].coordinates().y, Eq(q_earth));

  Length const q_probe = (trajectory.last().degrees_of_freedom().position() -
                         ICRFJ2000Equator::origin).coordinates().y;
  Speed const v_probe =
      trajectory.last().degrees_of_freedom().velocity().coordinates().x;
  std::vector<Displacement<ICRFJ2000Equator>> probe_positions;
  for (DiscreteTrajectory<ICRFJ2000Equator>::Iterator it =
           trajectory.Begin();
       it != trajectory.End();
       ++it) {
    probe_positions.push_back(it.degrees_of_freedom().position() -
                              ICRFJ2000Equator::origin);
  }
  EXPECT_THAT(probe_positions.size(), Eq(476));
  EXPECT_THAT(probe_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe, 259));
  EXPECT_THAT(probe_positions.back().coordinates().y,
              Eq(q_probe));

  Instant const old_t_max = ephemeris.t_max();
  EXPECT_THAT(trajectory.last().time(), Lt(old_t_max));
  EXPECT_FALSE(
      ephemeris.FlowWithAdaptiveStep(
          &trajectory,
          intrinsic_acceleration,
          t0_ + std::numeric_limits<double>::infinity() * Second,
          Ephemeris<ICRFJ2000Equator>::AdaptiveStepParameters(
              DormandElMikkawyPrince1986RKN434FM<Position<ICRFJ2000Equator>>(),
              kMaxSteps,
              1E-9 * Metre,
              2.6E-15 * Metre / Second),
          /*max_ephemeris_steps=*/0));
  EXPECT_THAT(ephemeris.t_max(), Eq(old_t_max));
  EXPECT_THAT(trajectory.last().time(), Eq(old_t_max));
}

// The Earth and two massless probes, similar to the previous test but flowing
// with a fixed step.
TEST_F(EphemerisTest, EarthTwoProbes) {
  Length const kDistance1 = 1E9 * Metre;
  Length const kDistance2 = 3E9 * Metre;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> initial_state;
  Position<ICRFJ2000Equator> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  bodies.erase(bodies.begin() + 1);
  initial_state.erase(initial_state.begin() + 1);

  MassiveBody const* const earth = bodies[0].get();
  Position<ICRFJ2000Equator> const earth_position =
      initial_state[0].position();
  Velocity<ICRFJ2000Equator> const earth_velocity =
      initial_state[0].velocity();

  Ephemeris<ICRFJ2000Equator>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              period / 100));

  MasslessBody probe1;
  DiscreteTrajectory<ICRFJ2000Equator> trajectory1;
  trajectory1.Append(t0_,
                     DegreesOfFreedom<ICRFJ2000Equator>(
                         earth_position + Vector<Length, ICRFJ2000Equator>(
                             {0 * Metre, kDistance1, 0 * Metre}),
                         earth_velocity));
  auto const intrinsic_acceleration1 =
      [earth, kDistance1](Instant const& t) {
        return Vector<Acceleration, ICRFJ2000Equator>(
            {0 * SIUnit<Acceleration>(),
             earth->gravitational_parameter() / (kDistance1 * kDistance1),
             0 * SIUnit<Acceleration>()});
      };

  MasslessBody probe2;
  DiscreteTrajectory<ICRFJ2000Equator> trajectory2;
  trajectory2.Append(t0_,
                     DegreesOfFreedom<ICRFJ2000Equator>(
                         earth_position + Vector<Length, ICRFJ2000Equator>(
                             {0 * Metre, -kDistance2, 0 * Metre}),
                         earth_velocity));
  auto const intrinsic_acceleration2 =
      [earth, kDistance2](Instant const& t) {
        return Vector<Acceleration, ICRFJ2000Equator>(
            {0 * SIUnit<Acceleration>(),
             -earth->gravitational_parameter() / (kDistance2 * kDistance2),
             0 * SIUnit<Acceleration>()});
      };

  ephemeris.FlowWithFixedStep(
      {&trajectory1, &trajectory2},
      {intrinsic_acceleration1, intrinsic_acceleration2},
      t0_ + period,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
          McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
          period / 1000));

  ContinuousTrajectory<ICRFJ2000Equator> const& earth_trajectory =
      *ephemeris.trajectory(earth);

  ContinuousTrajectory<ICRFJ2000Equator>::Hint hint;
  DegreesOfFreedom<ICRFJ2000Equator> const earth_degrees_of_freedom =
      earth_trajectory.EvaluateDegreesOfFreedom(t0_ + period, &hint);
  Length const q_earth = (earth_degrees_of_freedom.position() -
                          ICRFJ2000Equator::origin).coordinates().y;
  Speed const v_earth = earth_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<ICRFJ2000Equator>> earth_positions;
  for (int i = 0; i <= 100; ++i) {
    earth_positions.push_back(
        earth_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
            ICRFJ2000Equator::origin);
  }

  EXPECT_THAT(earth_positions.size(), Eq(101));
  EXPECT_THAT(earth_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v_earth, 551));
  EXPECT_THAT(earth_positions[25].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v_earth, 230));
  EXPECT_THAT(earth_positions[50].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v_earth, 413));
  EXPECT_THAT(earth_positions[75].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v_earth, 622));
  EXPECT_THAT(earth_positions[100].coordinates().y, Eq(q_earth));

  Length const q_probe1 = (trajectory1.last().degrees_of_freedom().position() -
                     ICRFJ2000Equator::origin).coordinates().y;
  Length const q_probe2 = (trajectory2.last().degrees_of_freedom().position() -
                     ICRFJ2000Equator::origin).coordinates().y;
  Speed const v_probe1 =
      trajectory1.last().degrees_of_freedom().velocity().coordinates().x;
  Speed const v_probe2 =
      trajectory2.last().degrees_of_freedom().velocity().coordinates().x;
  std::vector<Displacement<ICRFJ2000Equator>> probe1_positions;
  std::vector<Displacement<ICRFJ2000Equator>> probe2_positions;
  for (DiscreteTrajectory<ICRFJ2000Equator>::Iterator it =
           trajectory1.Begin();
       it != trajectory1.End();
       ++it) {
    probe1_positions.push_back(it.degrees_of_freedom().position() -
                               ICRFJ2000Equator::origin);
  }
  for (DiscreteTrajectory<ICRFJ2000Equator>::Iterator it =
           trajectory2.Begin();
       it != trajectory2.End();
       ++it) {
    probe2_positions.push_back(it.degrees_of_freedom().position() -
                               ICRFJ2000Equator::origin);
  }
#if defined(WE_LOVE_228)
  EXPECT_THAT(probe1_positions.size(), Eq(2));
  EXPECT_THAT(probe2_positions.size(), Eq(2));
#else
  EXPECT_THAT(probe1_positions.size(), Eq(1001));
  EXPECT_THAT(probe2_positions.size(), Eq(1001));
#endif
  EXPECT_THAT(probe1_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe1, 15));
  EXPECT_THAT(probe2_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe2, 0));
  EXPECT_THAT(probe1_positions.back().coordinates().y,
              Eq(q_probe1));
  EXPECT_THAT(probe2_positions.back().coordinates().y,
              Eq(q_probe2));
}

TEST_F(EphemerisTest, Спутник1ToСпутник2) {
  auto const at_спутник_1_launch =
      SolarSystemFactory::AtСпутник1Launch(
          SolarSystemFactory::Accuracy::kAllBodiesAndOblateness);
  auto const at_спутник_2_launch =
      SolarSystemFactory::AtСпутник2Launch(
          SolarSystemFactory::Accuracy::kAllBodiesAndOblateness);

  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(
          /*fitting_tolerance=*/5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              /*step=*/45 * Minute));

  ephemeris->Prolong(at_спутник_2_launch->epoch());

  // Upper bounds, tight to the nearest order of magnitude.
  static std::map<SolarSystemFactory::Index,
                  Angle> const expected_angle_error = {{}};
  static std::map<SolarSystemFactory::Index,
                  double> const expected_parent_distance_error = {{}};
  static std::map<SolarSystemFactory::Index,
                  double> const expected_parent_offset_error = {
      {SolarSystemFactory::kAriel, 1E-3},
      {SolarSystemFactory::kDione, 1E-3},
      {SolarSystemFactory::kIo, 1E-3},
      {SolarSystemFactory::kTethys, 1E-3},
      {SolarSystemFactory::kTitania, 1E-3},
      {SolarSystemFactory::kUmbriel, 1E-3},
      {SolarSystemFactory::kCharon, 1E-4},
      {SolarSystemFactory::kEuropa, 1E-4},
      {SolarSystemFactory::kOberon, 1E-4},
      {SolarSystemFactory::kRhea, 1E-4},
      {SolarSystemFactory::kTitan, 1E-4},
      {SolarSystemFactory::kTriton, 1E-4},
      {SolarSystemFactory::kGanymede, 1E-5},
      {SolarSystemFactory::kIapetus, 1E-5},
      {SolarSystemFactory::kMoon, 1E-5},  // What is this?
      {SolarSystemFactory::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystemFactory::kCallisto, 1E-7},
      {SolarSystemFactory::kVenus, 1E-7},
      {SolarSystemFactory::kEarth, 1E-8},
      {SolarSystemFactory::kSaturn, 1E-8},
      {SolarSystemFactory::kUranus, 1E-8},
      {SolarSystemFactory::kMars, 1E-9},
      {SolarSystemFactory::kPluto, 1E-9},
      {SolarSystemFactory::kJupiter, 1E-10},
      {SolarSystemFactory::kEris, 1E-12},
      {SolarSystemFactory::kNeptune, 1E-12}};
  static std::map<SolarSystemFactory::Index,
                  double> const expected_position_error = {
      {SolarSystemFactory::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystemFactory::kIo, 1E-6},
      {SolarSystemFactory::kTethys, 1E-6},
      {SolarSystemFactory::kDione, 1E-7},
      {SolarSystemFactory::kEuropa, 1E-7},
      {SolarSystemFactory::kMoon, 1E-7},
      {SolarSystemFactory::kOberon, 1E-7},
      {SolarSystemFactory::kRhea, 1E-7},
      {SolarSystemFactory::kTitan, 1E-7},
      {SolarSystemFactory::kTitania, 1E-7},
      {SolarSystemFactory::kUmbriel, 1E-7},
      {SolarSystemFactory::kVenus, 1E-7},
      {SolarSystemFactory::kAriel, 1E-8},
      {SolarSystemFactory::kEarth, 1E-8},
      {SolarSystemFactory::kGanymede, 1E-8},
      {SolarSystemFactory::kIapetus, 1E-8},
      {SolarSystemFactory::kSaturn, 1E-8},
      {SolarSystemFactory::kSun, 1E-8},
      {SolarSystemFactory::kTriton, 1E-8},
      {SolarSystemFactory::kUranus, 1E-8},
      {SolarSystemFactory::kCallisto, 1E-9},
      {SolarSystemFactory::kCharon, 1E-9},
      {SolarSystemFactory::kJupiter, 1E-9},
      {SolarSystemFactory::kMars, 1E-9},
      {SolarSystemFactory::kPluto, 1E-9},
      {SolarSystemFactory::kNeptune, 1E-12},
      {SolarSystemFactory::kEris, 1E-13}};
  static std::map<SolarSystemFactory::Index,
                  double> const expected_velocity_error = {
      {SolarSystemFactory::kDione, 1E-3},
      {SolarSystemFactory::kIo, 1E-3},
      {SolarSystemFactory::kTethys, 1E-3},
      {SolarSystemFactory::kAriel, 1E-4},
      {SolarSystemFactory::kEuropa, 1E-4},
      {SolarSystemFactory::kOberon, 1E-4},
      {SolarSystemFactory::kRhea, 1E-4},
      {SolarSystemFactory::kTitania, 1E-4},
      {SolarSystemFactory::kTriton, 1E-4},
      {SolarSystemFactory::kUmbriel, 1E-4},
      {SolarSystemFactory::kCharon, 1E-5},
      {SolarSystemFactory::kTitan, 1E-5},
      {SolarSystemFactory::kUranus, 1E-5},
      {SolarSystemFactory::kGanymede, 1E-6},
      {SolarSystemFactory::kIapetus, 1E-6},
      {SolarSystemFactory::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystemFactory::kMoon, 1E-6},
      {SolarSystemFactory::kPluto, 1E-6},
      {SolarSystemFactory::kSaturn, 1E-6},
      {SolarSystemFactory::kCallisto, 1E-7},
      {SolarSystemFactory::kEarth, 1E-7},
      {SolarSystemFactory::kJupiter, 1E-7},
      {SolarSystemFactory::kSun, 1E-7},
      {SolarSystemFactory::kVenus, 1E-7},
      {SolarSystemFactory::kMars, 1E-8},
      {SolarSystemFactory::kNeptune, 1E-8},
      {SolarSystemFactory::kEris, 1E-10}};

  for (int i = SolarSystemFactory::kSun;
       i <= SolarSystemFactory::kLastBody;
       ++i) {
    SolarSystemFactory::Index const index =
        static_cast<SolarSystemFactory::Index>(i);
    ContinuousTrajectory<ICRFJ2000Equator> const& trajectory =
        at_спутник_1_launch->trajectory(*ephemeris,
                                        SolarSystemFactory::name(i));
    DegreesOfFreedom<ICRFJ2000Equator> final_state =
        at_спутник_2_launch->initial_state(SolarSystemFactory::name(i));
    double const position_error = RelativeError(
        final_state.position() - kSolarSystemBarycentreEquator,
        trajectory.EvaluatePosition(at_спутник_2_launch->epoch(), nullptr) -
            kSolarSystemBarycentreEquator);
    double const velocity_error = RelativeError(
        final_state.velocity(),
        trajectory.EvaluateVelocity(at_спутник_2_launch->epoch(), nullptr));
    EXPECT_THAT(position_error, Lt(expected_position_error.at(index)))
        << SolarSystemFactory::name(i);
    EXPECT_THAT(position_error, Gt(expected_position_error.at(index) / 10.0))
        << SolarSystemFactory::name(i);
    EXPECT_THAT(velocity_error, Lt(expected_velocity_error.at(index)))
        << SolarSystemFactory::name(i);
    EXPECT_THAT(velocity_error, Gt(expected_velocity_error.at(index) / 10.0))
        << SolarSystemFactory::name(i);
    if (i != SolarSystemFactory::kSun) {
      // Look at the error in the position relative to the parent.
      Vector<Length, ICRFJ2000Equator> expected =
          final_state.position() -
          at_спутник_2_launch->initial_state(
              SolarSystemFactory::name(SolarSystemFactory::parent(i))).
              position();
      Vector<Length, ICRFJ2000Equator> actual =
          trajectory.EvaluatePosition(at_спутник_2_launch->epoch(), nullptr) -
          at_спутник_1_launch->trajectory(
              *ephemeris,
              SolarSystemFactory::name(SolarSystemFactory::parent(i))).
                  EvaluatePosition(at_спутник_2_launch->epoch(), nullptr);
      if (expected_angle_error.find(index) != expected_angle_error.end()) {
        Area const product_of_norms = expected.Norm() * actual.Norm();
        Angle const angle = ArcTan(
            Wedge(expected, actual).Norm() / product_of_norms,
            InnerProduct(expected, actual) / product_of_norms);
        EXPECT_THAT(angle / Degree,
                    Gt(expected_angle_error.at(index) / Degree * 0.9))
            << SolarSystemFactory::name(i);
        EXPECT_THAT(angle / Degree,
                    Lt(expected_angle_error.at(index) / Degree * 1.1))
            << SolarSystemFactory::name(i);
      }
      if (expected_parent_distance_error.find(index) !=
          expected_parent_distance_error.end()) {
        double const parent_distance_error = RelativeError(expected.Norm(),
                                                  actual.Norm());
        EXPECT_THAT(parent_distance_error,
                    Lt(expected_parent_distance_error.at(index)))
            << SolarSystemFactory::name(i);
        EXPECT_THAT(parent_distance_error,
                    Gt(expected_parent_distance_error.at(index) / 10.0))
            << SolarSystemFactory::name(i);
      }
      if (expected_parent_offset_error.find(index) !=
          expected_parent_offset_error.end()) {
        double const parent_offset_error =  RelativeError(expected, actual);
        EXPECT_THAT(parent_offset_error,
                    Lt(expected_parent_offset_error.at(index)))
            << SolarSystemFactory::name(i);
        EXPECT_THAT(parent_offset_error,
                    Gt(expected_parent_offset_error.at(index) / 10.0))
            << SolarSystemFactory::name(i);
      }
    }
  }
}

TEST_F(EphemerisTest, Serialization) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> initial_state;
  Position<ICRFJ2000Equator> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<ICRFJ2000Equator>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              period / 100));
  ephemeris.Prolong(t0_ + period);

  EXPECT_EQ(0, ephemeris.serialization_index_for_body(earth));
  EXPECT_EQ(1, ephemeris.serialization_index_for_body(moon));
  EXPECT_EQ(earth, ephemeris.body_for_serialization_index(0));
  EXPECT_EQ(moon, ephemeris.body_for_serialization_index(1));

  serialization::Ephemeris message;
  ephemeris.WriteToMessage(&message);

  auto const ephemeris_read =
      Ephemeris<ICRFJ2000Equator>::ReadFromMessage(message);
  MassiveBody const* const earth_read = ephemeris_read->bodies()[0];
  MassiveBody const* const moon_read = ephemeris_read->bodies()[1];

  EXPECT_EQ(0, ephemeris_read->serialization_index_for_body(earth_read));
  EXPECT_EQ(1, ephemeris_read->serialization_index_for_body(moon_read));
  EXPECT_EQ(earth_read, ephemeris_read->body_for_serialization_index(0));
  EXPECT_EQ(moon_read, ephemeris_read->body_for_serialization_index(1));

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
  Time const kDuration = 1 * Second;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> initial_state;

  auto earth = SolarSystem<ICRFJ2000Equator>::MakeMassiveBody(
      solar_system_.gravity_model_message("Earth"));
  Velocity<ICRFJ2000Equator> const v;
  Position<ICRFJ2000Equator> const q = ICRFJ2000Equator::origin;

  bodies.push_back(std::move(earth));
  initial_state.emplace_back(q, v);

  Position<ICRFJ2000Equator> const earth_position =
      initial_state[0].position();
  Velocity<ICRFJ2000Equator> const earth_velocity =
      initial_state[0].velocity();

  Ephemeris<ICRFJ2000Equator>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              kDuration / 100));

  MasslessBody elephant;
  DiscreteTrajectory<ICRFJ2000Equator> trajectory;
  trajectory.Append(t0_,
                    DegreesOfFreedom<ICRFJ2000Equator>(
                        earth_position + Vector<Length, ICRFJ2000Equator>(
                            {0 * Metre, 0 * Metre, kEarthPolarRadius}),
                        earth_velocity));

  ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<ICRFJ2000Equator>::kNoIntrinsicAcceleration,
      t0_ + kDuration,
      Ephemeris<ICRFJ2000Equator>::AdaptiveStepParameters(
          DormandElMikkawyPrince1986RKN434FM<Position<ICRFJ2000Equator>>(),
          kMaxSteps,
          1E-9 * Metre,
          2.6E-15 * Metre / Second),
          Ephemeris<ICRFJ2000Equator>::unlimited_max_ephemeris_steps);

  Speed const v_elephant_x =
      trajectory.last().degrees_of_freedom().velocity().coordinates().x;
  Speed const v_elephant_y =
      trajectory.last().degrees_of_freedom().velocity().coordinates().y;
  std::vector<Displacement<ICRFJ2000Equator>> elephant_positions;
  std::vector<Vector<Acceleration, ICRFJ2000Equator>> elephant_accelerations;
  for (DiscreteTrajectory<ICRFJ2000Equator>::Iterator it =
           trajectory.Begin();
       it != trajectory.End();
       ++it) {
    elephant_positions.push_back(it.degrees_of_freedom().position() -
                                 ICRFJ2000Equator::origin);
    elephant_accelerations.push_back(
        ephemeris.ComputeGravitationalAccelerationOnMasslessBody(
            &trajectory, it.time()));
  }

  // The small residual in x comes from the fact that the cosine of the
  // declination (90 degrees) is not exactly zero, so the axis of our Earth is
  // slightly tilted.  This greatly annoys the elephant.
  EXPECT_THAT(elephant_positions.size(), Eq(9));
  EXPECT_THAT(elephant_positions.back().coordinates().x,
              VanishesBefore(1 * Metre, 0));
  EXPECT_THAT(elephant_positions.back().coordinates().y,
              AlmostEquals(kDuration * v_elephant_y, 0));
  EXPECT_LT(RelativeError(elephant_positions.back().coordinates().z,
                          kEarthPolarRadius), 8E-7);

  EXPECT_THAT(elephant_accelerations.size(), Eq(9));
  EXPECT_THAT(elephant_accelerations.back().coordinates().x,
              VanishesBefore(1 * Metre / Second / Second, 0));
  EXPECT_THAT(elephant_accelerations.back().coordinates().y,
              AlmostEquals(0 * SIUnit<Acceleration>(), 0));
  EXPECT_LT(RelativeError(elephant_accelerations.back().coordinates().z,
                          -9.832 * SIUnit<Acceleration>()), 6.7E-6);
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
                                            1 * Metre,
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
  Position<World> const q0 = World::origin +
      Vector<Length, World>({0 * AstronomicalUnit,
                             0 * AstronomicalUnit,
                             0 * AstronomicalUnit});
  Position<World> const q1 = World::origin +
      Vector<Length, World>({1 * AstronomicalUnit,
                             0 * AstronomicalUnit,
                             0 * AstronomicalUnit});
  Position<World> const q2 = World::origin +
      Vector<Length, World>({1 * AstronomicalUnit,
                             0 * AstronomicalUnit,
                             1 * AstronomicalUnit});
  Position<World> const q3 = World::origin +
      Vector<Length, World>({0 * AstronomicalUnit,
                             0 * AstronomicalUnit,
                             1 * AstronomicalUnit});
  initial_state.emplace_back(q0, v);
  initial_state.emplace_back(q1, v);
  initial_state.emplace_back(q2, v);
  initial_state.emplace_back(q3, v);

  Ephemeris<World>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          5 * Milli(Metre),
          Ephemeris<World>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<World>>(),
              kDuration / 100));
  ephemeris.Prolong(t0_ + kDuration);

  Vector<Acceleration, World> actual_acceleration0 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b0, t0_);
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
              AlmostEquals(expected_acceleration0, 0, 6));

  Vector<Acceleration, World> actual_acceleration1 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b1, t0_);
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
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b2, t0_);
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
              AlmostEquals(expected_acceleration2, 0, 3));

  Vector<Acceleration, World> actual_acceleration3 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b3, t0_);
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

TEST_F(EphemerisTest, ComputeApsides) {
  Instant const t0;
  GravitationalParameter const μ = GravitationalConstant * SolarMass;
  auto const b = new MassiveBody(μ);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<World>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b));
  initial_state.emplace_back(World::origin, Velocity<World>());

  Ephemeris<World>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0,
          5 * Milli(Metre),
          Ephemeris<World>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<World>>(),
              1 * Hour));

  Displacement<World> r(
      {1 * AstronomicalUnit, 2 * AstronomicalUnit, 3 * AstronomicalUnit});
  Length const r_norm = r.Norm();
  Velocity<World> v({4 * Kilo(Metre) / Second,
                     5 * Kilo(Metre) / Second,
                     6 * Kilo(Metre) / Second});
  Speed const v_norm = v.Norm();

  Time const T = 2 * π * Sqrt(-(Pow<3>(r_norm) * Pow<2>(μ) /
                                Pow<3>(r_norm * Pow<2>(v_norm) - 2 * μ)));
  Length const a = -r_norm * μ / (r_norm * Pow<2>(v_norm) - 2 * μ);

  DiscreteTrajectory<World> trajectory;
  trajectory.Append(t0, DegreesOfFreedom<World>(World::origin + r, v));

  ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<World>::kNoIntrinsicAcceleration,
      t0 + 10 * JulianYear,
      Ephemeris<World>::AdaptiveStepParameters(
          DormandElMikkawyPrince1986RKN434FM<Position<World>>(),
          std::numeric_limits<std::int64_t>::max(),
          1E-3 * Metre,
          1E-3 * Metre / Second));

  auto const apsides =
      ephemeris.ComputeApsides(b, trajectory.Begin(), trajectory.End());
  std::experimental::optional<Instant> previous_time;
  std::experimental::optional<Position<World>> previous_position;
  int apsides_count = 0;
  for (auto it = apsides->Begin(); it != apsides->End(); ++it) {
    Instant const time = it.time();
    Position<World> const position = it.degrees_of_freedom().position();
    if (previous_time) {
      EXPECT_THAT(time - *previous_time,
                  AlmostEquals(0.5 * T, 103, 2143));
      EXPECT_THAT((position - *previous_position).Norm(),
                  AlmostEquals(2.0 * a, 0, 176));
    }
    previous_time = time;
    previous_position = position;
    ++apsides_count;
  }
  EXPECT_EQ(6, apsides_count);
}

}  // namespace physics
}  // namespace principia
