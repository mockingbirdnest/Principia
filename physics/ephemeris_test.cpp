
#include "physics/ephemeris.hpp"

#include <limits>
#include <map>
#include <optional>
#include <set>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/macros.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/constants.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {
namespace internal_ephemeris {

using astronomy::ICRS;
using base::not_null;
using geometry::Barycentre;
using geometry::AngularVelocity;
using geometry::Displacement;
using geometry::Frame;
using geometry::Rotation;
using geometry::Velocity;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using integrators::methods::McLachlanAtela1992Order4Optimal;
using integrators::methods::McLachlanAtela1992Order5Optimal;
using integrators::methods::Quinlan1999Order8A;
using quantities::Abs;
using quantities::ArcTan;
using quantities::Area;
using quantities::Mass;
using quantities::Pow;
using quantities::SIUnit;
using quantities::Sqrt;
using quantities::astronomy::JulianYear;
using quantities::astronomy::LunarDistance;
using quantities::astronomy::SolarMass;
using quantities::constants::GravitationalConstant;
using quantities::si::AstronomicalUnit;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::AbsoluteError;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystemFactory;
using testing_utilities::StatusIs;
using ::testing::AnyOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Ref;

namespace {

Length constexpr earth_polar_radius = 6356.8 * Kilo(Metre);
int constexpr max_steps = 1e6;
char constexpr big_name[] = "Big";
char constexpr small_name[] = "Small";

}  // namespace

class EphemerisTest : public testing::TestWithParam<FixedStepSizeIntegrator<
                          Ephemeris<ICRS>::NewtonianMotionEquation> const*> {
 protected:
  EphemerisTest()
      : solar_system_(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2433282_500000000.proto.txt"),
        t0_(solar_system_.epoch()) {}

  FixedStepSizeIntegrator<Ephemeris<ICRS>::NewtonianMotionEquation> const&
  integrator() {
    return *GetParam();
  }

  void SetUpEarthMoonSystem(
      std::vector<not_null<std::unique_ptr<MassiveBody const>>>& bodies,
      std::vector<DegreesOfFreedom<ICRS>>& initial_state,
      Position<ICRS>& centre_of_mass,
      Time& period) {
    // Make the bodies non-oblate so that the system can be computed explicitly.
    solar_system_.RemoveOblateness("Earth");
    solar_system_.RemoveOblateness("Moon");
    serialization::GravityModel::Body const earth_gravity_model =
        solar_system_.gravity_model_message("Earth");
    serialization::GravityModel::Body const moon_gravity_model =
        solar_system_.gravity_model_message("Moon");

    // Create the Moon before the Earth to exercise a bug caused by the order of
    // pointers differing from the order of bodies (don't ask).
    auto moon = SolarSystem<ICRS>::MakeMassiveBody(moon_gravity_model);
    auto earth = SolarSystem<ICRS>::MakeMassiveBody(earth_gravity_model);

    // The Earth-Moon system, roughly, with a circular orbit with velocities
    // in the centre-of-mass frame.
    Position<ICRS> const q1 =
        ICRS::origin + Displacement<ICRS>({0 * Metre, 0 * Metre, 0 * Metre});
    Position<ICRS> const q2 =
        ICRS::origin + Displacement<ICRS>({0 * Metre, 4e8 * Metre, 0 * Metre});
    Length const semi_major_axis = (q1 - q2).Norm();
    period = 2 * π *
             Sqrt(Pow<3>(semi_major_axis) / (earth->gravitational_parameter() +
                                             moon->gravitational_parameter()));
    centre_of_mass = Barycentre<Position<ICRS>, Mass>(
        {q1, q2}, {earth->mass(), moon->mass()});
    Velocity<ICRS> const v1({-2 * π * (q1 - centre_of_mass).Norm() / period,
                             0 * SIUnit<Speed>(),
                             0 * SIUnit<Speed>()});
    Velocity<ICRS> const v2({2 * π * (q2 - centre_of_mass).Norm() / period,
                             0 * SIUnit<Speed>(),
                             0 * SIUnit<Speed>()});

    bodies.push_back(std::move(earth));
    bodies.push_back(std::move(moon));
    initial_state.emplace_back(q1, v1);
    initial_state.emplace_back(q2, v2);
  }

  SolarSystem<ICRS> solar_system_;
  Instant t0_;
};

TEST_P(EphemerisTest, ProlongSpecialCases) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  Position<ICRS> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(bodies, initial_state, centre_of_mass, period);

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));
  EXPECT_THAT(ephemeris.planetary_integrator(), Ref(integrator()));

  EXPECT_EQ(astronomy::InfinitePast, ephemeris.t_max());
  EXPECT_EQ(astronomy::InfiniteFuture, ephemeris.t_min());

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

TEST_P(EphemerisTest, FlowWithAdaptiveStepSpecialCase) {
  Length const distance = 1e9 * Metre;
  Speed const velocity = 1e3 * Metre / Second;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  Position<ICRS> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(bodies, initial_state, centre_of_mass, period);

  Position<ICRS> const earth_position = initial_state[0].position();

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  MasslessBody probe;
  DiscreteTrajectory<ICRS> trajectory;
  trajectory.Append(
      t0_,
      DegreesOfFreedom<ICRS>(
          earth_position + Displacement<ICRS>({0 * Metre, distance, 0 * Metre}),
          Velocity<ICRS>({velocity, velocity, velocity})));

  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<ICRS>::NoIntrinsicAcceleration,
      t0_ + period,
      Ephemeris<ICRS>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Position<ICRS>>(),
          max_steps,
          1e-9 * Metre,
          2.6e-15 * Metre / Second),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps,
      /*last_point_only=*/false));
  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<ICRS>::NoIntrinsicAcceleration,
      trajectory.last().time(),
      Ephemeris<ICRS>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Position<ICRS>>(),
          max_steps,
          1e-9 * Metre,
          2.6e-15 * Metre / Second),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps,
      /*last_point_only=*/false));
}

// The canonical Earth-Moon system, tuned to produce circular orbits.
TEST_P(EphemerisTest, EarthMoon) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  Position<ICRS> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(bodies, initial_state, centre_of_mass, period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  ephemeris.Prolong(t0_ + period);

  ContinuousTrajectory<ICRS> const& earth_trajectory =
      *ephemeris.trajectory(earth);
  ContinuousTrajectory<ICRS> const& moon_trajectory =
      *ephemeris.trajectory(moon);

  std::vector<Displacement<ICRS>> earth_positions;
  for (int i = 0; i <= 100; ++i) {
    earth_positions.push_back(
      earth_trajectory.EvaluatePosition(t0_ + i * period / 100) -
          centre_of_mass);
  }
  EXPECT_THAT(earth_positions.size(), Eq(101));
  EXPECT_THAT(Abs(earth_positions[25].coordinates().y), Lt(6e-4 * Metre));
  EXPECT_THAT(Abs(earth_positions[50].coordinates().x), Lt(7e-3 * Metre));
  EXPECT_THAT(Abs(earth_positions[75].coordinates().y), Lt(2e-2 * Metre));
  EXPECT_THAT(Abs(earth_positions[100].coordinates().x), Lt(3e-2 * Metre));

  std::vector<Displacement<ICRS>> moon_positions;
  for (int i = 0; i <= 100; ++i) {
    moon_positions.push_back(
      moon_trajectory.EvaluatePosition(t0_ + i * period / 100) -
          centre_of_mass);
  }
  EXPECT_THAT(moon_positions.size(), Eq(101));
  EXPECT_THAT(Abs(moon_positions[25].coordinates().y), Lt(5e-2 * Metre));
  EXPECT_THAT(Abs(moon_positions[50].coordinates().x), Lt(6e-1 * Metre));
  EXPECT_THAT(Abs(moon_positions[75].coordinates().y), Lt(2 * Metre));
  EXPECT_THAT(Abs(moon_positions[100].coordinates().x), Lt(2 * Metre));
}

// Test the behavior of ForgetBefore on the Earth-Moon system.
TEST_P(EphemerisTest, ForgetBefore) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  Position<ICRS> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(bodies, initial_state, centre_of_mass, period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 10));

  ephemeris.Prolong(t0_ + 16 * period);

  ContinuousTrajectory<ICRS> const& earth_trajectory =
      *ephemeris.trajectory(earth);
  ContinuousTrajectory<ICRS> const& moon_trajectory =
      *ephemeris.trajectory(moon);

  Instant t_max = ephemeris.t_max();
  EXPECT_EQ(t0_ + 16 * period, t_max);
  EXPECT_EQ(t_max, earth_trajectory.t_max());
  EXPECT_EQ(t_max, moon_trajectory.t_max());

  ephemeris.ForgetBefore(t0_ + 3 * period);
  EXPECT_EQ(t0_ + 3 * period, ephemeris.t_min());
  EXPECT_EQ(t0_ + 3 * period, earth_trajectory.t_min());
  EXPECT_EQ(t0_ + 3 * period, moon_trajectory.t_min());
}

// The Moon alone.  It moves in straight line.
TEST_P(EphemerisTest, Moon) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  Position<ICRS> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(bodies, initial_state, centre_of_mass, period);

  bodies.erase(bodies.begin());
  initial_state.erase(initial_state.begin());

  MassiveBody const* const moon = bodies[0].get();

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  ephemeris.Prolong(t0_ + period);

  ContinuousTrajectory<ICRS> const& moon_trajectory =
      *ephemeris.trajectory(moon);

  DegreesOfFreedom<ICRS> const moon_degrees_of_freedom =
      moon_trajectory.EvaluateDegreesOfFreedom(t0_ + period);
  Length const q =
      (moon_degrees_of_freedom.position() - ICRS::origin).coordinates().y;
  Speed const v = moon_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<ICRS>> moon_positions;
  for (int i = 0; i <= 100; ++i) {
    moon_positions.push_back(
        moon_trajectory.EvaluatePosition(t0_ + i * period / 100) -
            ICRS::origin);
  }

  EXPECT_THAT(moon_positions.size(), Eq(101));
  EXPECT_THAT(moon_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v, 343, 345));
  EXPECT_THAT(moon_positions[25].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v, 151, 153));
  EXPECT_THAT(moon_positions[50].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v, 512, 516));
  EXPECT_THAT(moon_positions[75].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v, 401, 403));
  EXPECT_THAT(moon_positions[100].coordinates().y, Eq(q));
}

// The Earth and a massless probe 1 billion meters away, with the same velocity,
// and an acceleration which exactly compensates gravitational attraction.  Both
// bodies move in straight lines.
TEST_P(EphemerisTest, EarthProbe) {
  Length const distance = 1e9 * Metre;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  Position<ICRS> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(bodies, initial_state, centre_of_mass, period);

  bodies.erase(bodies.begin() + 1);
  initial_state.erase(initial_state.begin() + 1);

  MassiveBody const* const earth = bodies[0].get();
  Position<ICRS> const earth_position = initial_state[0].position();
  Velocity<ICRS> const earth_velocity = initial_state[0].velocity();

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  MasslessBody probe;
  DiscreteTrajectory<ICRS> trajectory;
  trajectory.Append(t0_,
                    DegreesOfFreedom<ICRS>(
                        earth_position + Vector<Length, ICRS>(
                                             {0 * Metre, distance, 0 * Metre}),
                        earth_velocity));
  auto const intrinsic_acceleration = [earth, distance](Instant const& t) {
    return Vector<Acceleration, ICRS>(
        {0 * SIUnit<Acceleration>(),
         earth->gravitational_parameter() / (distance * distance),
         0 * SIUnit<Acceleration>()});
  };

  ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      intrinsic_acceleration,
      t0_ + period,
      Ephemeris<ICRS>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Position<ICRS>>(),
          max_steps,
          1e-9 * Metre,
          2.6e-15 * Metre / Second),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps,
      /*last_point_only=*/false);

  ContinuousTrajectory<ICRS> const& earth_trajectory =
      *ephemeris.trajectory(earth);

  DegreesOfFreedom<ICRS> const earth_degrees_of_freedom =
      earth_trajectory.EvaluateDegreesOfFreedom(t0_ + period);
  Length const q_earth =
      (earth_degrees_of_freedom.position() - ICRS::origin).coordinates().y;
  Speed const v_earth = earth_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<ICRS>> earth_positions;
  for (int i = 0; i <= 100; ++i) {
    earth_positions.push_back(
        earth_trajectory.EvaluatePosition(t0_ + i * period / 100) -
            ICRS::origin);
  }

  EXPECT_THAT(earth_positions.size(), Eq(101));
  EXPECT_THAT(earth_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v_earth, 539, 541));
  EXPECT_THAT(earth_positions[25].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v_earth, 241, 242));
  EXPECT_THAT(earth_positions[50].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v_earth, 403, 405));
  EXPECT_THAT(earth_positions[75].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v_earth, 633, 635));
  EXPECT_THAT(earth_positions[100].coordinates().y, Eq(q_earth));

  Length const q_probe = (trajectory.last().degrees_of_freedom().position() -
                          ICRS::origin).coordinates().y;
  Speed const v_probe =
      trajectory.last().degrees_of_freedom().velocity().coordinates().x;
  std::vector<Displacement<ICRS>> probe_positions;
  for (DiscreteTrajectory<ICRS>::Iterator it = trajectory.Begin();
       it != trajectory.End();
       ++it) {
    probe_positions.push_back(it.degrees_of_freedom().position() -
                              ICRS::origin);
  }
  // The solution is a line, so the rounding errors dominate.  Different
  // libms result in different errors and thus different numbers of steps.
  EXPECT_THAT(probe_positions.size(),
              AnyOf(Eq(410),    // MSVC Release
                    Eq(421),    // MSVC Debug
                    Eq(446)));  // Clang Linux
  EXPECT_THAT(probe_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe, 222, 259));
  EXPECT_THAT(probe_positions.back().coordinates().y,
              Eq(q_probe));

  Instant const old_t_max = ephemeris.t_max();
  EXPECT_THAT(trajectory.last().time(), Lt(old_t_max));
  EXPECT_THAT(ephemeris.FlowWithAdaptiveStep(
                  &trajectory,
                  intrinsic_acceleration,
                  t0_ + std::numeric_limits<double>::infinity() * Second,
                  Ephemeris<ICRS>::AdaptiveStepParameters(
                      EmbeddedExplicitRungeKuttaNyströmIntegrator<
                          DormandالمكاوىPrince1986RKN434FM,
                          Position<ICRS>>(),
                      max_steps,
                      1e-9 * Metre,
                      2.6e-15 * Metre / Second),
                  /*max_ephemeris_steps=*/0,
                  /*last_point_only=*/false),
              StatusIs(Error::DEADLINE_EXCEEDED));
  EXPECT_THAT(ephemeris.t_max(), Eq(old_t_max));
  EXPECT_THAT(trajectory.last().time(), Eq(old_t_max));
}

// The Earth and two massless probes, similar to the previous test but flowing
// with a fixed step.
TEST_P(EphemerisTest, EarthTwoProbes) {
  Length const distance_1 = 1e9 * Metre;
  Length const distance_2 = 3e9 * Metre;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  Position<ICRS> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(bodies, initial_state, centre_of_mass, period);

  bodies.erase(bodies.begin() + 1);
  initial_state.erase(initial_state.begin() + 1);

  MassiveBody const* const earth = bodies[0].get();
  Position<ICRS> const earth_position = initial_state[0].position();
  Velocity<ICRS> const earth_velocity = initial_state[0].velocity();

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  MasslessBody probe1;
  DiscreteTrajectory<ICRS> trajectory1;
  trajectory1.Append(t0_,
                     DegreesOfFreedom<ICRS>(
                         earth_position + Vector<Length, ICRS>(
                             {0 * Metre, distance_1, 0 * Metre}),
                         earth_velocity));
  auto const intrinsic_acceleration1 = [earth, distance_1](Instant const& t) {
    return Vector<Acceleration, ICRS>(
        {0 * SIUnit<Acceleration>(),
         earth->gravitational_parameter() / (distance_1 * distance_1),
         0 * SIUnit<Acceleration>()});
  };

  MasslessBody probe2;
  DiscreteTrajectory<ICRS> trajectory2;
  trajectory2.Append(t0_,
                     DegreesOfFreedom<ICRS>(
                         earth_position + Vector<Length, ICRS>(
                             {0 * Metre, -distance_2, 0 * Metre}),
                         earth_velocity));
  auto const intrinsic_acceleration2 = [earth, distance_2](Instant const& t) {
    return Vector<Acceleration, ICRS>(
        {0 * SIUnit<Acceleration>(),
         -earth->gravitational_parameter() / (distance_2 * distance_2),
         0 * SIUnit<Acceleration>()});
  };

  auto instance = ephemeris.NewInstance(
      {&trajectory1, &trajectory2},
      {intrinsic_acceleration1, intrinsic_acceleration2},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 1000));
  EXPECT_OK(ephemeris.FlowWithFixedStep(t0_ + period, *instance));

  ContinuousTrajectory<ICRS> const& earth_trajectory =
      *ephemeris.trajectory(earth);

  DegreesOfFreedom<ICRS> const earth_degrees_of_freedom =
      earth_trajectory.EvaluateDegreesOfFreedom(t0_ + period);
  Length const q_earth =
      (earth_degrees_of_freedom.position() - ICRS::origin).coordinates().y;
  Speed const v_earth = earth_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<ICRS>> earth_positions;
  for (int i = 0; i <= 100; ++i) {
    earth_positions.push_back(
        earth_trajectory.EvaluatePosition(t0_ + i * period / 100) -
            ICRS::origin);
  }

  EXPECT_THAT(earth_positions.size(), Eq(101));
  EXPECT_THAT(earth_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v_earth, 539, 541));
  EXPECT_THAT(earth_positions[25].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v_earth, 241, 242));
  EXPECT_THAT(earth_positions[50].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v_earth, 403, 405));
  EXPECT_THAT(earth_positions[75].coordinates().y, Eq(q_earth));
  EXPECT_THAT(earth_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v_earth, 633, 635));
  EXPECT_THAT(earth_positions[100].coordinates().y, Eq(q_earth));

  Length const q_probe1 = (trajectory1.last().degrees_of_freedom().position() -
                     ICRS::origin).coordinates().y;
  Length const q_probe2 = (trajectory2.last().degrees_of_freedom().position() -
                     ICRS::origin).coordinates().y;
  Speed const v_probe1 =
      trajectory1.last().degrees_of_freedom().velocity().coordinates().x;
  Speed const v_probe2 =
      trajectory2.last().degrees_of_freedom().velocity().coordinates().x;
  std::vector<Displacement<ICRS>> probe1_positions;
  std::vector<Displacement<ICRS>> probe2_positions;
  for (DiscreteTrajectory<ICRS>::Iterator it = trajectory1.Begin();
       it != trajectory1.End();
       ++it) {
    probe1_positions.push_back(it.degrees_of_freedom().position() -
                               ICRS::origin);
  }
  for (DiscreteTrajectory<ICRS>::Iterator it = trajectory2.Begin();
       it != trajectory2.End();
       ++it) {
    probe2_positions.push_back(it.degrees_of_freedom().position() -
                               ICRS::origin);
  }
  EXPECT_THAT(probe1_positions.size(), Eq(1001));
  EXPECT_THAT(probe2_positions.size(), Eq(1001));
  EXPECT_THAT(probe1_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe1, 25, 70));
  EXPECT_THAT(probe2_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe2, 2, 13));
  EXPECT_THAT(probe1_positions.back().coordinates().y,
              Eq(q_probe1));
  EXPECT_THAT(probe2_positions.back().coordinates().y,
              Eq(q_probe2));
}

TEST_P(EphemerisTest, Serialization) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  Position<ICRS> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(bodies, initial_state, centre_of_mass, period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));
  ephemeris.Prolong(t0_ + period);

  EXPECT_EQ(0, ephemeris.serialization_index_for_body(earth));
  EXPECT_EQ(1, ephemeris.serialization_index_for_body(moon));
  EXPECT_EQ(earth, ephemeris.body_for_serialization_index(0));
  EXPECT_EQ(moon, ephemeris.body_for_serialization_index(1));

  serialization::Ephemeris message;
  ephemeris.WriteToMessage(&message);

  auto const ephemeris_read = Ephemeris<ICRS>::ReadFromMessage(message);
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
    EXPECT_EQ(
        ephemeris.trajectory(earth)->EvaluateDegreesOfFreedom(time),
        ephemeris_read->trajectory(earth_read)->EvaluateDegreesOfFreedom(time));
    EXPECT_EQ(
        ephemeris.trajectory(moon)->EvaluateDegreesOfFreedom(time),
        ephemeris_read->trajectory(moon_read)->EvaluateDegreesOfFreedom(time));
  }

  serialization::Ephemeris second_message;
  ephemeris_read->WriteToMessage(&second_message);
  EXPECT_THAT(message, EqualsProto(second_message));
}

// The gravitational acceleration on an elephant located at the pole.
TEST_P(EphemerisTest, ComputeGravitationalAccelerationMasslessBody) {
  Time const duration = 1 * Second;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;

  auto earth = SolarSystem<ICRS>::MakeMassiveBody(
      solar_system_.gravity_model_message("Earth"));
  Velocity<ICRS> const v;
  Position<ICRS> const q = ICRS::origin;

  bodies.push_back(std::move(earth));
  initial_state.emplace_back(q, v);

  Position<ICRS> const earth_position = initial_state[0].position();
  Velocity<ICRS> const earth_velocity = initial_state[0].velocity();

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), duration / 100));

  // The elephant's initial position and velocity.
  DiscreteTrajectory<ICRS> trajectory;
  trajectory.Append(t0_,
                    DegreesOfFreedom<ICRS>(
                        earth_position + Vector<Length, ICRS>(
                            {0 * Metre, 0 * Metre, earth_polar_radius}),
                        earth_velocity));

  ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<ICRS>::NoIntrinsicAcceleration,
      t0_ + duration,
      Ephemeris<ICRS>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Position<ICRS>>(),
          max_steps,
          1e-9 * Metre,
          2.6e-15 * Metre / Second),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps,
      /*last_point_only=*/false);

  Speed const v_elephant_y =
      trajectory.last().degrees_of_freedom().velocity().coordinates().y;
  std::vector<Displacement<ICRS>> elephant_positions;
  std::vector<Vector<Acceleration, ICRS>> elephant_accelerations;
  for (DiscreteTrajectory<ICRS>::Iterator it = trajectory.Begin();
       it != trajectory.End();
       ++it) {
    elephant_positions.push_back(it.degrees_of_freedom().position() -
                                 ICRS::origin);
    elephant_accelerations.push_back(
        ephemeris.ComputeGravitationalAccelerationOnMasslessBody(
            &trajectory, it.time()));
  }

  // The small residual in x comes from the fact that the cosine of the
  // declination (90 degrees) is not exactly zero, so the axis of our Earth is
  // slightly tilted.  Also, the geopotential is not rotationally symmetrical,
  // so there is a tiny residual in y too.  This greatly annoys the elephant.
  EXPECT_THAT(elephant_positions.size(), Eq(9));
  EXPECT_THAT(elephant_positions.back().coordinates().x,
              IsNear(-9.8e-19 * Metre));
  EXPECT_THAT(elephant_positions.back().coordinates().y,
              AnyOf(IsNear(-2.3e-31 * Metre), Eq(0 * Metre)));
  EXPECT_LT(RelativeError(elephant_positions.back().coordinates().z,
                          earth_polar_radius), 8e-7);

  EXPECT_THAT(elephant_accelerations.size(), Eq(9));
  EXPECT_THAT(elephant_accelerations.back().coordinates().x,
              IsNear(-2.0e-18 * Metre / Second / Second));
  EXPECT_THAT(elephant_accelerations.back().coordinates().y,
              AnyOf(IsNear(-2.7e-30 * Metre / Second / Second),
                    Eq(0 * Metre / Second / Second)));
  EXPECT_LT(RelativeError(elephant_accelerations.back().coordinates().z,
                          -9.832 * SIUnit<Acceleration>()), 6.7e-6);
}

// An apple located a bit above the pole collides with the ground.
TEST_P(EphemerisTest, CollisionDetection) {
  Time const short_duration = 1 * Second;
  // A naïve computation would have 1.428 s.
  Time const long_duration = 1.431 * Second;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;

  auto earth = SolarSystem<ICRS>::MakeMassiveBody(
      solar_system_.gravity_model_message("Earth"));
  Length const& earth_mean_radius = earth->mean_radius();
  Velocity<ICRS> const v;
  Position<ICRS> const q = ICRS::origin;

  bodies.push_back(std::move(earth));
  initial_state.emplace_back(q, v);

  Position<ICRS> const earth_position = initial_state[0].position();
  Velocity<ICRS> const earth_velocity = initial_state[0].velocity();

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), short_duration / 100));

  // The apple's initial position and velocity
  DiscreteTrajectory<ICRS> trajectory;
  trajectory.Append(
      t0_,
      DegreesOfFreedom<ICRS>(
          earth_position +
              Vector<Length, ICRS>(
                  {0 * Metre, 0 * Metre, earth_mean_radius + 10 * Metre}),
          earth_velocity));
  auto const instance = ephemeris.NewInstance(
      {&trajectory},
      Ephemeris<ICRS>::NoIntrinsicAccelerations,
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                             Position<ICRS>>(),
          1e-3 * Second));

  EXPECT_OK(ephemeris.FlowWithFixedStep(t0_ + short_duration, *instance));
  EXPECT_THAT(ephemeris.FlowWithFixedStep(t0_ + long_duration, *instance),
              StatusIs(Error::OUT_OF_RANGE));
}

TEST_P(EphemerisTest, ComputeGravitationalAccelerationMassiveBody) {
  Time const duration = 1 * Second;
  double const j2 = 1e6;
  Length const radius = 1 * LunarDistance;

  Mass const m0 = 1 * SolarMass;
  Mass const m1 = 2 * SolarMass;
  Mass const m2 = 3 * SolarMass;
  Mass const m3 = 4 * SolarMass;

  auto const b0 = new OblateBody<ICRS>(
      m0,
      RotatingBody<ICRS>::Parameters(1 * Metre,
                                                 1 * Radian,
                                                 t0_,
                                                 4 * Radian / Second,
                                                 0 * Radian,
                                                 π / 2 * Radian),
      OblateBody<ICRS>::Parameters(j2, radius));
  auto const b1 = new MassiveBody(m1);
  auto const b2 = new MassiveBody(m2);
  auto const b3 = new MassiveBody(m3);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b0));
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b1));
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b2));
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b3));

  Velocity<ICRS> const v({0 * SIUnit<Speed>(),
                                      0 * SIUnit<Speed>(),
                                      0 * SIUnit<Speed>()});
  Position<ICRS> const q0 = ICRS::origin +
      Vector<Length, ICRS>({0 * AstronomicalUnit,
                                        0 * AstronomicalUnit,
                                        0 * AstronomicalUnit});
  Position<ICRS> const q1 = ICRS::origin +
      Vector<Length, ICRS>({1 * AstronomicalUnit,
                                        0 * AstronomicalUnit,
                                        0 * AstronomicalUnit});
  Position<ICRS> const q2 = ICRS::origin +
      Vector<Length, ICRS>({1 * AstronomicalUnit,
                                        0 * AstronomicalUnit,
                                        1 * AstronomicalUnit});
  Position<ICRS> const q3 = ICRS::origin +
      Vector<Length, ICRS>({0 * AstronomicalUnit,
                                        0 * AstronomicalUnit,
                                        1 * AstronomicalUnit});
  initial_state.emplace_back(q0, v);
  initial_state.emplace_back(q1, v);
  initial_state.emplace_back(q2, v);
  initial_state.emplace_back(q3, v);

  Ephemeris<ICRS> ephemeris(
      std::move(bodies),
      initial_state,
      t0_,
      5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(integrator(), duration / 100));
  ephemeris.Prolong(t0_ + duration);

  Vector<Acceleration, ICRS> actual_acceleration0 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b0, t0_);
  Vector<Acceleration, ICRS> expected_acceleration0 =
      GravitationalConstant * (m1 * (q1 - q0) / Pow<3>((q1 - q0).Norm()) +
                               m2 * (q2 - q0) / Pow<3>((q2 - q0).Norm()) +
                               m3 * (q3 - q0) / Pow<3>((q3 - q0).Norm())) +
      Vector<Acceleration, ICRS>(
          {(1.5 * m1 - (9 / Sqrt(512)) * m2) * GravitationalConstant *
               Pow<2>(radius) * j2 / Pow<4>((q0 - q1).Norm()),
           0 * SIUnit<Acceleration>(),
           (-3 * m3 + (3 / Sqrt(512)) * m2) * GravitationalConstant *
               Pow<2>(radius) * j2 / Pow<4>((q0 - q1).Norm())});
  EXPECT_THAT(actual_acceleration0,
              AlmostEquals(expected_acceleration0, 0, 6));

  Vector<Acceleration, ICRS> actual_acceleration1 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b1, t0_);
  Vector<Acceleration, ICRS> expected_acceleration1 =
      GravitationalConstant * (m0 * (q0 - q1) / Pow<3>((q0 - q1).Norm()) +
                               m2 * (q2 - q1) / Pow<3>((q2 - q1).Norm()) +
                               m3 * (q3 - q1) / Pow<3>((q3 - q1).Norm())) +
      Vector<Acceleration, ICRS>(
          {-1.5 * GravitationalConstant * m0 * Pow<2>(radius) * j2 /
               Pow<4>((q0 - q1).Norm()),
           0 * SIUnit<Acceleration>(),
           0 * SIUnit<Acceleration>()});
  EXPECT_THAT(actual_acceleration1,
              AlmostEquals(expected_acceleration1, 0, 4));

  Vector<Acceleration, ICRS> actual_acceleration2 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b2, t0_);
  Vector<Acceleration, ICRS> expected_acceleration2 =
      GravitationalConstant * (m0 * (q0 - q2) / Pow<3>((q0 - q2).Norm()) +
                               m1 * (q1 - q2) / Pow<3>((q1 - q2).Norm()) +
                               m3 * (q3 - q2) / Pow<3>((q3 - q2).Norm())) +
      Vector<Acceleration, ICRS>(
          {(9 / Sqrt(512)) * GravitationalConstant * m0 *
               Pow<2>(radius) * j2 / Pow<4>((q0 - q1).Norm()),
           0 * SIUnit<Acceleration>(),
           (-3 / Sqrt(512)) * GravitationalConstant * m0 *
               Pow<2>(radius) * j2 / Pow<4>((q0 - q1).Norm())});
  EXPECT_THAT(actual_acceleration2,
              AlmostEquals(expected_acceleration2, 0, 3));

  Vector<Acceleration, ICRS> actual_acceleration3 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b3, t0_);
  Vector<Acceleration, ICRS> expected_acceleration3 =
      GravitationalConstant * (m0 * (q0 - q3) / Pow<3>((q0 - q3).Norm()) +
                               m1 * (q1 - q3) / Pow<3>((q1 - q3).Norm()) +
                               m2 * (q2 - q3) / Pow<3>((q2 - q3).Norm())) +
      Vector<Acceleration, ICRS>(
          {0 * SIUnit<Acceleration>(),
           0 * SIUnit<Acceleration>(),
           3 * GravitationalConstant * m0 * Pow<2>(radius) * j2 /
               Pow<4>((q0 - q1).Norm())});
  EXPECT_THAT(actual_acceleration3,
              AlmostEquals(expected_acceleration3, 0, 4));
}

TEST_P(EphemerisTest, ComputeApsidesContinuousTrajectory) {
  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "test_gravity_model_two_bodies.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "test_initial_state_two_bodies_elliptical.proto.txt");

  Length const fitting_tolerance = 1 * Milli(Metre);
  Instant const t0 = solar_system.epoch();
  Time const T =
      16000 * π / (Sqrt(7) * std::pow(73 - 8 * Sqrt(35), 1.5)) * Second;
  Length const a = 400 / (73 - 8 * Sqrt(35)) * Kilo(Metre);
  double const e = (7 + 8 * Sqrt(35)) / 80;
  Speed v_apoapsis = (-1631 + 348 * Sqrt(35)) / 1460 * Kilo(Metre) / Second;
  Speed v_periapsis = Sqrt(7 * (87 + 8 * Sqrt(35))) / 20 * Kilo(Metre) / Second;

  auto ephemeris = solar_system.MakeEphemeris(
      fitting_tolerance,
      Ephemeris<ICRS>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<McLachlanAtela1992Order4Optimal,
                                                Position<ICRS>>(),
          /*step=*/10 * Milli(Second)));
  ephemeris->Prolong(t0 + 10 * T);

  MassiveBody const* const big =
      solar_system.massive_body(*ephemeris, big_name);
  MassiveBody const* const small =
      solar_system.massive_body(*ephemeris, small_name);
  DiscreteTrajectory<ICRS> apoapsides1;
  DiscreteTrajectory<ICRS> apoapsides2;
  DiscreteTrajectory<ICRS> periapsides1;
  DiscreteTrajectory<ICRS> periapsides2;
  ephemeris->ComputeApsides(big,
                            small,
                            apoapsides1,
                            periapsides1,
                            apoapsides2,
                            periapsides2);

  EXPECT_EQ(apoapsides1.Size(), apoapsides2.Size());
  EXPECT_EQ(periapsides1.Size(), periapsides2.Size());

  EXPECT_EQ(10, apoapsides1.Size());
  EXPECT_EQ(10, periapsides1.Size());

  std::optional<Instant> previous_time;
  std::set<Instant> all_times;
  for (auto it1 = apoapsides1.Begin(), it2 = apoapsides2.Begin();
       it1 != apoapsides1.End() && it2 != apoapsides2.End();
       ++it1, ++it2) {
    Instant const time = it1.time();
    all_times.emplace(time);
    Displacement<ICRS> const displacement =
        it1.degrees_of_freedom().position() -
        it2.degrees_of_freedom().position();
    EXPECT_LT(AbsoluteError(displacement.Norm(), (1 + e) * a),
              1.9e-5 * fitting_tolerance);
    if (previous_time) {
      EXPECT_LT(AbsoluteError(time - *previous_time, T),
                0.11 * fitting_tolerance / v_apoapsis);
    }
    previous_time = time;
  }

  previous_time = std::nullopt;
  for (auto it1 = periapsides1.Begin(), it2 = periapsides2.Begin();
       it1 != periapsides1.End() && it2 != periapsides2.End();
       ++it1, ++it2) {
    Instant const time = it1.time();
    all_times.emplace(time);
    Displacement<ICRS> const displacement =
        it1.degrees_of_freedom().position() -
        it2.degrees_of_freedom().position();
    EXPECT_LT(AbsoluteError(displacement.Norm(), (1 - e) * a),
              5.3e-3 * fitting_tolerance);
    if (previous_time) {
      EXPECT_LT(AbsoluteError(time - *previous_time, T),
                2.1 * fitting_tolerance / v_periapsis);
    }
    previous_time = time;
  }

  previous_time = std::nullopt;
  for (Instant const& time : all_times) {
    if (previous_time) {
      EXPECT_LT(AbsoluteError(time - *previous_time, 0.5 * T),
                2.3 * fitting_tolerance / (v_apoapsis + v_periapsis));
    }
    previous_time = time;
  }
}

INSTANTIATE_TEST_CASE_P(
    AllEphemerisTests,
    EphemerisTest,
    ::testing::Values(
        &SymplecticRungeKuttaNyströmIntegrator<McLachlanAtela1992Order5Optimal,
                                               Position<ICRS>>(),
        &SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                            Position<ICRS>>()));

}  // namespace internal_ephemeris
}  // namespace physics
}  // namespace principia
