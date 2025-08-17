#include "physics/ephemeris.hpp"

#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/space.hpp"
#include "gipfeli/gipfeli.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/integrators.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/elementary_functions.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rotating_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/arithmetic.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {

using ::testing::AllOf;
using ::testing::AnyOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Ref;
using namespace principia::astronomy::_frames;
using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_generalized_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::integrators::_symplectic_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::numerics::_elementary_functions;
using namespace principia::physics::_continuous_trajectory;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_oblate_body;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics;
using namespace principia::testing_utilities::_numerics_matchers;
using namespace principia::testing_utilities::_vanishes_before;
using namespace std::chrono_literals;

namespace {

int constexpr max_steps = 1'000'000;
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
    solar_system_.LimitOblatenessToDegree("Earth", /*max_degree=*/0);
    solar_system_.LimitOblatenessToDegree("Moon", /*max_degree=*/0);
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
    centre_of_mass = Barycentre({q1, q2}, {earth->mass(), moon->mass()});
    Velocity<ICRS> const v1({-2 * π * (q1 - centre_of_mass).Norm() / period,
                             0 * si::Unit<Speed>,
                             0 * si::Unit<Speed>});
    Velocity<ICRS> const v2({2 * π * (q2 - centre_of_mass).Norm() / period,
                             0 * si::Unit<Speed>,
                             0 * si::Unit<Speed>});

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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));
  EXPECT_THAT(ephemeris.planetary_integrator(), Ref(integrator()));

  EXPECT_EQ(InfinitePast, ephemeris.t_max());
  EXPECT_EQ(InfiniteFuture, ephemeris.t_min());

  EXPECT_TRUE(ephemeris.empty());
  EXPECT_OK(ephemeris.Prolong(t0_ + period));
  EXPECT_FALSE(ephemeris.empty());
  EXPECT_EQ(t0_, ephemeris.t_min());
  EXPECT_LE(t0_ + period, ephemeris.t_max());
  Instant const t_max = ephemeris.t_max();

  EXPECT_OK(ephemeris.Prolong(t0_ + period / 2));
  EXPECT_EQ(t_max, ephemeris.t_max());

  Instant const last_t = Barycentre({t0_ + period, t_max});
  ephemeris.Prolong(last_t).IgnoreError();
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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  MasslessBody probe;
  DiscreteTrajectory<ICRS> trajectory;
  EXPECT_OK(trajectory.Append(
      t0_,
      DegreesOfFreedom<ICRS>(
          earth_position + Displacement<ICRS>({0 * Metre, distance, 0 * Metre}),
          Velocity<ICRS>({velocity, velocity, velocity}))));

  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<ICRS>::NoIntrinsicAcceleration,
      t0_ + period,
      Ephemeris<ICRS>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          max_steps,
          1e-9 * Metre,
          2.6e-15 * Metre / Second),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps));
  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<ICRS>::NoIntrinsicAcceleration,
      trajectory.back().time,
      Ephemeris<ICRS>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          max_steps,
          1e-9 * Metre,
          2.6e-15 * Metre / Second),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps));
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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  EXPECT_OK(ephemeris.Prolong(t0_ + period));

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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  EXPECT_OK(ephemeris.Prolong(t0_ + period));

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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  MasslessBody probe;
  DiscreteTrajectory<ICRS> trajectory;
  EXPECT_OK(trajectory.Append(
      t0_,
      DegreesOfFreedom<ICRS>(
          earth_position +
              Vector<Length, ICRS>({0 * Metre, distance, 0 * Metre}),
          earth_velocity)));
  auto const intrinsic_acceleration = [earth, distance](Instant const& t) {
    return Vector<Acceleration, ICRS>(
        {0 * si::Unit<Acceleration>,
         earth->gravitational_parameter() / (distance * distance),
         0 * si::Unit<Acceleration>});
  };

  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      intrinsic_acceleration,
      t0_ + period,
      Ephemeris<ICRS>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          max_steps,
          1e-9 * Metre,
          2.6e-15 * Metre / Second),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps));

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

  Length const q_probe = (trajectory.back().degrees_of_freedom.position() -
                          ICRS::origin).coordinates().y;
  Speed const v_probe =
      trajectory.back().degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<ICRS>> probe_positions;
  for (auto const& [time, degrees_of_freedom] : trajectory) {
    probe_positions.push_back(degrees_of_freedom.position() - ICRS::origin);
  }
  // The solution is a line, so the rounding errors dominate.  The size depends
  // on the test parameter (the fixed-step integrator of the ephemeris), on the
  // libm, including its use of FMA (via the std::pow in the adaptive step), and
  // on the use of FMA in polynomial evaluation.
  EXPECT_THAT(probe_positions.size(),
              AnyOf(Eq(358),    // MSVC no FMA/0
                    Eq(366),    // Clang Linux
                    Eq(373),    // MSVC all FMA/0
                    Eq(406),    // MSVC all FMA/1
                    Eq(420),    // MSVC FMA in libm only/1
                    Eq(421)));  // MSVC no FMA/1
  EXPECT_THAT(probe_positions.back().coordinates().x,
              AlmostEquals(1.00 * period * v_probe, 220, 266));
  EXPECT_THAT(probe_positions.back().coordinates().y,
              Eq(q_probe));

  Instant const old_t_max = ephemeris.t_max();
  EXPECT_THAT(trajectory.back().time, Lt(old_t_max));
  EXPECT_THAT(ephemeris.FlowWithAdaptiveStep(
                  &trajectory,
                  intrinsic_acceleration,
                  t0_ + std::numeric_limits<double>::infinity() * Second,
                  Ephemeris<ICRS>::AdaptiveStepParameters(
                      EmbeddedExplicitRungeKuttaNyströmIntegrator<
                          DormandالمكاوىPrince1986RKN434FM,
                          Ephemeris<ICRS>::NewtonianMotionEquation>(),
                      max_steps,
                      1e-9 * Metre,
                      2.6e-15 * Metre / Second),
                  /*max_ephemeris_steps=*/0),
              StatusIs(absl::StatusCode::kDeadlineExceeded));
  EXPECT_THAT(ephemeris.t_max(), Eq(old_t_max));
  EXPECT_THAT(trajectory.back().time, Eq(old_t_max));
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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));

  MasslessBody probe1;
  DiscreteTrajectory<ICRS> trajectory1;
  EXPECT_OK(trajectory1.Append(
      t0_,
      DegreesOfFreedom<ICRS>(
          earth_position +
              Vector<Length, ICRS>({0 * Metre, distance_1, 0 * Metre}),
          earth_velocity)));
  auto const intrinsic_acceleration1 = [earth, distance_1](Instant const& t) {
    return Vector<Acceleration, ICRS>(
        {0 * si::Unit<Acceleration>,
         earth->gravitational_parameter() / (distance_1 * distance_1),
         0 * si::Unit<Acceleration>});
  };

  MasslessBody probe2;
  DiscreteTrajectory<ICRS> trajectory2;
  EXPECT_OK(trajectory2.Append(
      t0_,
      DegreesOfFreedom<ICRS>(
          earth_position +
              Vector<Length, ICRS>({0 * Metre, -distance_2, 0 * Metre}),
          earth_velocity)));
  auto const intrinsic_acceleration2 = [earth, distance_2](Instant const& t) {
    return Vector<Acceleration, ICRS>(
        {0 * si::Unit<Acceleration>,
         -earth->gravitational_parameter() / (distance_2 * distance_2),
         0 * si::Unit<Acceleration>});
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

  Length const q_probe1 =
      (trajectory1.back().degrees_of_freedom.position() -
       ICRS::origin).coordinates().y;
  Length const q_probe2 =
      (trajectory2.back().degrees_of_freedom.position() -
       ICRS::origin).coordinates().y;
  Speed const v_probe1 =
      trajectory1.back().degrees_of_freedom.velocity().coordinates().x;
  Speed const v_probe2 =
      trajectory2.back().degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<ICRS>> probe1_positions;
  std::vector<Displacement<ICRS>> probe2_positions;
  for (auto const& [time, degrees_of_freedom] : trajectory1) {
    probe1_positions.push_back(degrees_of_freedom.position() - ICRS::origin);
  }
  for (auto const& [time, degrees_of_freedom] : trajectory2) {
    probe2_positions.push_back(degrees_of_freedom.position() - ICRS::origin);
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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), period / 100));
  EXPECT_OK(ephemeris.Prolong(t0_ + period));

  EXPECT_EQ(0, ephemeris.serialization_index_for_body(earth));
  EXPECT_EQ(1, ephemeris.serialization_index_for_body(moon));
  EXPECT_EQ(earth, ephemeris.body_for_serialization_index(0));
  EXPECT_EQ(moon, ephemeris.body_for_serialization_index(1));

  serialization::Ephemeris message;
  ephemeris.WriteToMessage(&message);

  auto const ephemeris_read = Ephemeris<ICRS>::ReadFromMessage(
      /*desired_t_min=*/InfiniteFuture,
      message);
  // After deserialization, the client must prolong as needed.
  EXPECT_OK(ephemeris_read->Prolong(ephemeris.t_max()));

  MassiveBody const* const earth_read = ephemeris_read->bodies()[0];
  MassiveBody const* const moon_read = ephemeris_read->bodies()[1];

  EXPECT_EQ(0, ephemeris_read->serialization_index_for_body(earth_read));
  EXPECT_EQ(1, ephemeris_read->serialization_index_for_body(moon_read));
  EXPECT_EQ(earth_read, ephemeris_read->body_for_serialization_index(0));
  EXPECT_EQ(moon_read, ephemeris_read->body_for_serialization_index(1));

  EXPECT_EQ(ephemeris.t_min(), ephemeris_read->t_min());
  for (Instant time = ephemeris.t_min();
       time <= ephemeris.t_max();
       time += (ephemeris.t_max() - ephemeris.t_min()) / 100) {
    EXPECT_OK(ephemeris_read->Prolong(time));
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

  solar_system_.LimitOblatenessToDegree("Earth", /*max_degree=*/2);
  solar_system_.LimitOblatenessToZonal("Earth");
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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), duration / 100));

  // The elephant's initial position and velocity.
  DiscreteTrajectory<ICRS> trajectory;
  EXPECT_OK(trajectory.Append(
      t0_,
      DegreesOfFreedom<ICRS>(
          earth_position + Vector<Length, ICRS>(
                               {0 * Metre, 0 * Metre, TerrestrialPolarRadius}),
          earth_velocity)));

  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<ICRS>::NoIntrinsicAcceleration,
      t0_ + duration,
      Ephemeris<ICRS>::GeneralizedAdaptiveStepParameters(
          EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
              Fine1987RKNG34,
              Ephemeris<ICRS>::GeneralizedNewtonianMotionEquation>(),
          max_steps,
          1e-9 * Metre,
          2.6e-15 * Metre / Second),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps));

  std::vector<Displacement<ICRS>> elephant_positions;
  std::vector<Vector<Acceleration, ICRS>> elephant_accelerations;
  for (auto const& [time, degrees_of_freedom] : trajectory) {
    elephant_positions.push_back(degrees_of_freedom.position() - ICRS::origin);
    elephant_accelerations.push_back(
        ephemeris.ComputeGravitationalAccelerationOnMasslessBody(
            &trajectory, time));
  }

  // The small residual in x comes from the fact that the cosine of the
  // declination (90 degrees) is not exactly zero, so the axis of our Earth is
  // slightly tilted.  Also, the geopotential is not rotationally symmetrical,
  // so there is a tiny residual in y too.  This greatly annoys the elephant.
  EXPECT_THAT(elephant_positions.size(), Eq(5));
  EXPECT_THAT(elephant_positions.back().coordinates().x,
              IsNear(-9.8e-19_(1) * Metre));
  EXPECT_THAT(elephant_positions.back().coordinates().y,
              AnyOf(IsNear(6.0e-35_(1) * Metre), Eq(0 * Metre)));
  EXPECT_LT(RelativeError(elephant_positions.back().coordinates().z,
                          TerrestrialPolarRadius), 8e-7);

  EXPECT_THAT(elephant_accelerations.size(), Eq(5));
  EXPECT_THAT(elephant_accelerations.back().coordinates().x,
              IsNear(-2.0e-18_(1) * Metre / Second / Second));
  EXPECT_THAT(elephant_accelerations.back().coordinates().y,
              AnyOf(IsNear(1.2e-34_(1) * Metre / Second / Second),
                    Eq(0 * Metre / Second / Second)));
  EXPECT_LT(RelativeError(elephant_accelerations.back().coordinates().z,
                          -9.832 * si::Unit<Acceleration>), 6.7e-6);
}

#if !defined(_DEBUG)
// An apple located a bit above the pole collides with the ground.
TEST_P(EphemerisTest, CollisionDetection) {
  Time const short_duration = 1 * Second;
  // A naïve computation would have 360.4 s.
  Time const long_duration = 354.8 * Second;
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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), short_duration / 100));

  // The apple's initial position and velocity.
  DiscreteTrajectory<ICRS> trajectory;
  EXPECT_OK(trajectory.Append(
      t0_,
      DegreesOfFreedom<ICRS>(
          earth_position +
              Vector<Length, ICRS>(
                  {0 * Metre, 0 * Metre, earth_mean_radius + 10 * Metre}),
          earth_velocity)));
  auto const instance =
      ephemeris.NewInstance({&trajectory},
                            Ephemeris<ICRS>::NoIntrinsicAccelerations,
                            Ephemeris<ICRS>::FixedStepParameters(
                                SymmetricLinearMultistepIntegrator<
                                    Quinlan1999Order8A,
                                    Ephemeris<ICRS>::NewtonianMotionEquation>(),
                                1e-3 * Second));

  EXPECT_OK(ephemeris.FlowWithFixedStep(t0_ + short_duration, *instance));
  EXPECT_THAT(ephemeris.FlowWithFixedStep(t0_ + long_duration, *instance),
              StatusIs(absl::StatusCode::kOutOfRange));
}
#endif

TEST_P(EphemerisTest, ComputeJacobianOnMassiveBody) {
  SolarSystem<ICRS> const solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  Instant const j2000 = solar_system_2000.epoch();

  // The most convenient way to compute the acceleration field near the centre
  // of the earth is to remove the earth.
  SolarSystem<ICRS> solar_system_2000_without_earth(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  solar_system_2000_without_earth.RemoveMassiveBody("Earth");

  Ephemeris<ICRS>::AccuracyParameters const accuracy_parameters {
    /*fitting_tolerance-*/ 1 * Milli(Metre),
    /*geopotential_tolerance=*/0x1p-24};
  Ephemeris<ICRS>::FixedStepParameters const fixed_step_parameters(
      SymplecticRungeKuttaNyströmIntegrator<
          McLachlanAtela1992Order4Optimal,
          Ephemeris<ICRS>::NewtonianMotionEquation>(),
      /*step=*/10 * Minute);

  auto ephemeris = solar_system_2000.MakeEphemeris(
      accuracy_parameters,
      fixed_step_parameters);
  CHECK_OK(ephemeris->Prolong(j2000));

  auto ephemeris_without_earth = solar_system_2000_without_earth.MakeEphemeris(
      accuracy_parameters,
      fixed_step_parameters);
  CHECK_OK(ephemeris_without_earth->Prolong(j2000));

  auto const earth = ephemeris->bodies()[solar_system_2000.index("Earth")];
  auto const earth_trajectory = ephemeris->trajectory(earth);
  Position<ICRS> const earth_position =
      earth_trajectory->EvaluatePosition(j2000);

  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> shift_distribution(1, 2);
  for (int i = 0; i < 1000; ++i) {
    Displacement<ICRS> const shift_x(
        {shift_distribution(random) * Metre,
         0 * Metre,
         0 * Metre});
    Displacement<ICRS> const shift_y(
        {0 * Metre,
         shift_distribution(random) * Metre,
         0 * Metre});
    Displacement<ICRS> const shift_z(
        {0 * Metre,
         0 * Metre,
         shift_distribution(random) * Metre});

    auto const acceleration_base =
        ephemeris_without_earth->ComputeGravitationalAccelerationOnMasslessBody(
            earth_position, j2000).coordinates();
    auto const acceleration_shift_x =
        ephemeris_without_earth->ComputeGravitationalAccelerationOnMasslessBody(
            earth_position + shift_x, j2000).coordinates();
    auto const acceleration_shift_y =
        ephemeris_without_earth->ComputeGravitationalAccelerationOnMasslessBody(
            earth_position + shift_y, j2000).coordinates();
    auto const acceleration_shift_z =
        ephemeris_without_earth->ComputeGravitationalAccelerationOnMasslessBody(
            earth_position + shift_z, j2000).coordinates();
    auto const finite_difference_jacobian_coordinates =
        R3x3Matrix<Inverse<Square<Time>>>(
            {(acceleration_shift_x.x - acceleration_base.x) /
                 shift_x.coordinates().x,
             (acceleration_shift_y.x - acceleration_base.x) /
                 shift_y.coordinates().y,
             (acceleration_shift_z.x - acceleration_base.x) /
                 shift_z.coordinates().z},
            {(acceleration_shift_x.y - acceleration_base.y) /
                 shift_x.coordinates().x,
             (acceleration_shift_y.y - acceleration_base.y) /
                 shift_y.coordinates().y,
             (acceleration_shift_z.y - acceleration_base.y) /
                 shift_z.coordinates().z},
            {(acceleration_shift_x.z - acceleration_base.z) /
                 shift_x.coordinates().x,
             (acceleration_shift_y.z - acceleration_base.z) /
                 shift_y.coordinates().y,
             (acceleration_shift_z.z - acceleration_base.z) /
                 shift_z.coordinates().z});
    auto const actual_jacobian =
        ephemeris->ComputeJacobianOnMassiveBody(earth, j2000);

    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        EXPECT_THAT(finite_difference_jacobian_coordinates(r, c),
                    RelativeErrorFrom(actual_jacobian.coordinates()(r, c),
                                      AllOf(Gt(2.3e-10), Lt(8.7e-5))))
            << r << " " << c;
      }
    }
  }
}

TEST_P(EphemerisTest, ComputeGravitationalJerkOnMasslessBody) {
  SolarSystem<ICRS> const solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  Instant const j2000 = solar_system_2000.epoch();

  auto ephemeris = solar_system_2000.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance-*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<
              McLachlanAtela1992Order4Optimal,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/10 * Minute));
  CHECK_OK(ephemeris->Prolong(j2000));

  auto const earth = ephemeris->bodies()[solar_system_2000.index("Earth")];
  auto const earth_trajectory = ephemeris->trajectory(earth);
  auto const earth_degrees_of_freedom =
      earth_trajectory->EvaluateDegreesOfFreedom(j2000);

  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> delay_distribution(1, 5);
  std::uniform_real_distribution<double> length_distribution(-1e9, 1e9);
  std::uniform_real_distribution<double> speed_distribution(-1e3, 1e3);
  for (int i = 0; i < 10; ++i) {
    Instant previous_t = j2000;
    DiscreteTrajectory<ICRS> vessel_trajectory;
    Displacement<ICRS> const displacement(
        {length_distribution(random) * Metre,
         length_distribution(random) * Metre,
         length_distribution(random) * Metre});
    Velocity<ICRS> const velocity(
        {speed_distribution(random) * Metre / Second,
         speed_distribution(random) * Metre / Second,
         speed_distribution(random) * Metre / Second});
    CHECK_OK(vessel_trajectory.Append(
        previous_t,
        earth_degrees_of_freedom +
            RelativeDegreesOfFreedom<ICRS>(displacement, velocity)));
    auto const instance = ephemeris->NewInstance(
        {&vessel_trajectory},
        Ephemeris<ICRS>::NoIntrinsicAccelerations,
        Ephemeris<ICRS>::FixedStepParameters(
            integrator(), 10 * Second));
    for (Instant t = j2000 + 1 * Minute; t < j2000 + 100 * Minute;
         t += delay_distribution(random) * Minute) {
      CHECK_OK(ephemeris->Prolong(t));
      CHECK_OK(ephemeris->FlowWithFixedStep(t + 10 * Minute, *instance));
      auto const previous_acceleration =
          ephemeris->ComputeGravitationalAccelerationOnMasslessBody(
              vessel_trajectory.EvaluatePosition(previous_t), previous_t);
      auto const acceleration =
          ephemeris->ComputeGravitationalAccelerationOnMasslessBody(
              vessel_trajectory.EvaluatePosition(t), t);
      auto const finite_difference_jerk =
          (acceleration - previous_acceleration) / (t - previous_t);
      auto const actual_jerk =
          ephemeris->ComputeGravitationalJerkOnMasslessBody(
              vessel_trajectory.EvaluateDegreesOfFreedom(t), t);

      EXPECT_THAT(
          finite_difference_jerk,
          Componentwise(RelativeErrorFrom(actual_jerk.coordinates().x,
                                          AllOf(Gt(9.6e-7), Lt(1.9e-3))),
                        RelativeErrorFrom(actual_jerk.coordinates().y,
                                          AllOf(Gt(6.0e-6), Lt(2.2e-3))),
                        RelativeErrorFrom(actual_jerk.coordinates().z,
                                          AllOf(Gt(1.8e-6), Lt(4.3e-3)))));

      previous_t = t;
    }
  }
}

TEST_P(EphemerisTest, ComputeGravitationalJerkOnMassiveBody) {
  SolarSystem<ICRS> const solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  Instant const j2000 = solar_system_2000.epoch();

  auto ephemeris = solar_system_2000.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance-*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<
              McLachlanAtela1992Order4Optimal,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/10 * Minute));
  CHECK_OK(ephemeris->Prolong(j2000));

  auto const earth = ephemeris->bodies()[solar_system_2000.index("Earth")];

  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> delay_distribution(1, 5);
  Instant previous_t = j2000;
  auto previous_acceleration =
      ephemeris->ComputeGravitationalAccelerationOnMassiveBody(earth,
                                                               previous_t);
  for (Instant t = j2000 + 1 * Minute;
       t < j2000 + 1000 * Minute;
       t += delay_distribution(random) * Minute) {
    CHECK_OK(ephemeris->Prolong(t));
    auto const acceleration =
        ephemeris->ComputeGravitationalAccelerationOnMassiveBody(earth, t);
    auto const finite_difference_jerk =
        (acceleration - previous_acceleration) / (t - previous_t);
    auto const actual_jerk =
        ephemeris->ComputeGravitationalJerkOnMassiveBody(earth, t);

    EXPECT_THAT(
        finite_difference_jerk,
        Componentwise(RelativeErrorFrom(actual_jerk.coordinates().x,
                                        AllOf(Gt(6.6e-7), Lt(4.7e-6))),
                      RelativeErrorFrom(actual_jerk.coordinates().y,
                                        AllOf(Gt(6.3e-5), Lt(3.4e-4))),
                      RelativeErrorFrom(actual_jerk.coordinates().z,
                                        AllOf(Gt(5.9e-5), Lt(3.1e-4)))));

    previous_t = t;
    previous_acceleration = acceleration;
  }
}

TEST_P(EphemerisTest, ComputeGravitationalAccelerationOnMassiveBody) {
  Time const duration = 1 * Second;
  double const j2 = 1e6;
  Length const radius = 1 * TerrestrialEquatorialRadius;

  auto const μ0 = 1 * SolarGravitationalParameter;
  auto const μ1 = 2 * SolarGravitationalParameter;
  auto const μ2 = 3 * SolarGravitationalParameter;
  auto const μ3 = 4 * SolarGravitationalParameter;

  auto const b0 =
      new OblateBody<ICRS>(μ0,
                           RotatingBody<ICRS>::Parameters(1 * Metre,
                                                          1 * Radian,
                                                          t0_,
                                                          4 * Radian / Second,
                                                          0 * Radian,
                                                          π / 2 * Radian),
                           OblateBody<ICRS>::Parameters(j2, radius));
  auto const b1 = new MassiveBody(μ1);
  auto const b2 = new MassiveBody(μ2);
  auto const b3 = new MassiveBody(μ3);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<ICRS>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b0));
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b1));
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b2));
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b3));

  Velocity<ICRS> const v({0 * si::Unit<Speed>,
                          0 * si::Unit<Speed>,
                          0 * si::Unit<Speed>});
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
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(integrator(), duration / 100));
  EXPECT_OK(ephemeris.Prolong(t0_ + duration));

  Vector<Acceleration, ICRS> actual_acceleration0 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b0, t0_);
  Vector<Acceleration, ICRS> expected_acceleration0 =
      (μ1 * (q1 - q0) / Pow<3>((q1 - q0).Norm()) +
       μ2 * (q2 - q0) / Pow<3>((q2 - q0).Norm()) +
       μ3 * (q3 - q0) / Pow<3>((q3 - q0).Norm())) +
      Vector<Acceleration, ICRS>(
          {(1.5 * μ1 - (9 / Sqrt(512)) * μ2) * Pow<2>(radius) * j2 /
               Pow<4>((q0 - q1).Norm()),
           0 * si::Unit<Acceleration>,
           (-3 * μ3 + (3 / Sqrt(512)) * μ2) * Pow<2>(radius) * j2 /
               Pow<4>((q0 - q1).Norm())});
  EXPECT_THAT(
      actual_acceleration0,
      Componentwise(AlmostEquals(expected_acceleration0.coordinates().x, 1),
                    VanishesBefore(expected_acceleration0.coordinates().x, 0),
                    AlmostEquals(expected_acceleration0.coordinates().z, 2)));

  Vector<Acceleration, ICRS> actual_acceleration1 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b1, t0_);
  Vector<Acceleration, ICRS> expected_acceleration1 =
      μ0 * (q0 - q1) / Pow<3>((q0 - q1).Norm()) +
      μ2 * (q2 - q1) / Pow<3>((q2 - q1).Norm()) +
      μ3 * (q3 - q1) / Pow<3>((q3 - q1).Norm()) +
      Vector<Acceleration, ICRS>(
          {-1.5 * μ0 * Pow<2>(radius) * j2 / Pow<4>((q0 - q1).Norm()),
           0 * si::Unit<Acceleration>,
           0 * si::Unit<Acceleration>});
  EXPECT_THAT(
      actual_acceleration1,
      Componentwise(AlmostEquals(expected_acceleration1.coordinates().x, 2),
                    VanishesBefore(expected_acceleration1.coordinates().x, 0),
                    AlmostEquals(expected_acceleration1.coordinates().z, 0)));

  Vector<Acceleration, ICRS> actual_acceleration2 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b2, t0_);
  Vector<Acceleration, ICRS> expected_acceleration2 =
      μ0 * (q0 - q2) / Pow<3>((q0 - q2).Norm()) +
      μ1 * (q1 - q2) / Pow<3>((q1 - q2).Norm()) +
      μ3 * (q3 - q2) / Pow<3>((q3 - q2).Norm()) +
      Vector<Acceleration, ICRS>({(9 / Sqrt(512)) * μ0 * Pow<2>(radius) * j2 /
                                      Pow<4>((q0 - q1).Norm()),
                                  0 * si::Unit<Acceleration>,
                                  (-3 / Sqrt(512)) * μ0 * Pow<2>(radius) * j2 /
                                      Pow<4>((q0 - q1).Norm())});
  EXPECT_THAT(
      actual_acceleration2,
      Componentwise(AlmostEquals(expected_acceleration2.coordinates().x, 4),
                    VanishesBefore(expected_acceleration2.coordinates().x, 0),
                    AlmostEquals(expected_acceleration2.coordinates().z, 1)));

  Vector<Acceleration, ICRS> actual_acceleration3 =
      ephemeris.ComputeGravitationalAccelerationOnMassiveBody(b3, t0_);
  Vector<Acceleration, ICRS> expected_acceleration3 =
      μ0 * (q0 - q3) / Pow<3>((q0 - q3).Norm()) +
      μ1 * (q1 - q3) / Pow<3>((q1 - q3).Norm()) +
      μ2 * (q2 - q3) / Pow<3>((q2 - q3).Norm()) +
      Vector<Acceleration, ICRS>(
          {0 * si::Unit<Acceleration>,
           0 * si::Unit<Acceleration>,
           3 * μ0 * Pow<2>(radius) * j2 / Pow<4>((q0 - q1).Norm())});
  EXPECT_THAT(
      actual_acceleration3,
      Componentwise(AlmostEquals(expected_acceleration3.coordinates().x, 3),
                    VanishesBefore(expected_acceleration3.coordinates().x, 0),
                    AlmostEquals(expected_acceleration3.coordinates().z, 0)));
}

TEST_P(EphemerisTest, ComputeGravitationalPotential) {
  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  Instant const j2000 = solar_system_2000.epoch();

  auto ephemeris = solar_system_2000.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance-*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<
              McLachlanAtela1992Order4Optimal,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/10 * Minute));
  CHECK_OK(ephemeris->Prolong(j2000));

  auto const earth_trajectory = ephemeris->trajectory(
      ephemeris->bodies()[solar_system_2000.index("Earth")]);
  Position<ICRS> const earth_position =
      earth_trajectory->EvaluatePosition(j2000);

  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> length_distribution(-1e7, 1e7);
  std::uniform_real_distribution<double> shift_distribution(1, 2);
  for (int i = 0; i < 1000; ++i) {
    Displacement<ICRS> const displacement(
        {length_distribution(random) * Metre,
         length_distribution(random) * Metre,
         length_distribution(random) * Metre});
    Displacement<ICRS> const shift_x(
        {shift_distribution(random) * Metre,
         0 * Metre,
         0 * Metre});
    Displacement<ICRS> const shift_y(
        {0 * Metre,
         shift_distribution(random) * Metre,
         0 * Metre});
    Displacement<ICRS> const shift_z(
        {0 * Metre,
         0 * Metre,
         shift_distribution(random) * Metre});

    SpecificEnergy const potential_base =
        ephemeris->ComputeGravitationalPotential(
            earth_position + displacement, j2000);
    SpecificEnergy const potential_shift_x =
        ephemeris->ComputeGravitationalPotential(
            earth_position + displacement + shift_x, j2000);
    SpecificEnergy const potential_shift_y =
        ephemeris->ComputeGravitationalPotential(
            earth_position + displacement + shift_y, j2000);
    SpecificEnergy const potential_shift_z =
        ephemeris->ComputeGravitationalPotential(
            earth_position + displacement + shift_z, j2000);
    auto const finite_difference_acceleration = Vector<Acceleration, ICRS>(
        {-(potential_shift_x - potential_base) /
             shift_x.coordinates().x,
         -(potential_shift_y - potential_base) /
             shift_y.coordinates().y,
         -(potential_shift_z - potential_base) /
             shift_z.coordinates().z});
    Vector<Acceleration, ICRS> const actual_acceleration =
        ephemeris->ComputeGravitationalAccelerationOnMasslessBody(
            earth_position + displacement, j2000);

    EXPECT_THAT(
        finite_difference_acceleration,
        RelativeErrorFrom(actual_acceleration, AllOf(Gt(1.1e-7), Lt(9.0e-6))));
  }
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
      /*accuracy_parameters=*/{fitting_tolerance,
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<
              McLachlanAtela1992Order4Optimal,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/10 * Milli(Second)));
  EXPECT_OK(ephemeris->Prolong(t0 + 10 * T));

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

  EXPECT_EQ(apoapsides1.size(), apoapsides2.size());
  EXPECT_EQ(periapsides1.size(), periapsides2.size());

  EXPECT_EQ(10, apoapsides1.size());
  EXPECT_EQ(10, periapsides1.size());

  std::optional<Instant> previous_time;
  std::set<Instant> all_times;
  for (auto it1 = apoapsides1.begin(), it2 = apoapsides2.begin();
       it1 != apoapsides1.end() && it2 != apoapsides2.end();
       ++it1, ++it2) {
    auto const& [time1, degrees_of_freedom1] = *it1;
    auto const& [time2, degrees_of_freedom2] = *it2;
    EXPECT_EQ(time1, time2);
    all_times.emplace(time1);
    Displacement<ICRS> const displacement =
        degrees_of_freedom1.position() - degrees_of_freedom2.position();
    EXPECT_LT(AbsoluteError(displacement.Norm(), (1 + e) * a),
              1.9e-5 * fitting_tolerance);
    if (previous_time) {
      EXPECT_LT(AbsoluteError(time1 - *previous_time, T),
                0.11 * fitting_tolerance / v_apoapsis);
    }
    previous_time = time1;
  }

  previous_time = std::nullopt;
  for (auto it1 = periapsides1.begin(), it2 = periapsides2.begin();
       it1 != periapsides1.end() && it2 != periapsides2.end();
       ++it1, ++it2) {
    auto const& [time1, degrees_of_freedom1] = *it1;
    auto const& [time2, degrees_of_freedom2] = *it2;
    EXPECT_EQ(time1, time2);
    all_times.emplace(time1);
    Displacement<ICRS> const displacement =
        degrees_of_freedom1.position() - degrees_of_freedom2.position();
    EXPECT_LT(AbsoluteError(displacement.Norm(), (1 - e) * a),
              5.3e-3 * fitting_tolerance);
    if (previous_time) {
      EXPECT_LT(AbsoluteError(time1 - *previous_time, T),
                2.1 * fitting_tolerance / v_periapsis);
    }
    previous_time = time1;
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

#if !defined(_DEBUG)
// This trajectory is similar to the second trajectory in the first save in
// #2400.  It exhibits oscillations with a period close to 5600 s and its
// downsampling period alternates between 120 and 130 s.
TEST(EphemerisTestNoFixture, DiscreteTrajectoryCompression) {
  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2433282_500000000.proto.txt");

  auto ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      /*fixed_step_parameters=*/{
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/10 * Minute});

  Instant const t0 = Instant() - 1.323698948906726e9 * Second;
  Instant const t1 = Instant() - 1.323595217786725e9 * Second;
  Position<ICRS> const q0 =
      Position<ICRS>{} + Displacement<ICRS>({-7.461169719467950e10 * Metre,
                                             1.165327644733623e11 * Metre,
                                             5.049935298178532e10 * Metre});
  Velocity<ICRS> const p0({-30169.49165384562 * Metre / Second,
                           -11880.03394238412 * Metre / Second,
                           200.2551546021678 * Metre / Second});

  MasslessBody probe;
  DiscreteTrajectory<ICRS> trajectory1;
  auto& segment1 = trajectory1.segments().front();
  segment1.SetDownsampling({.max_dense_intervals = 10'000,
                            .tolerance = 10 * Metre});
  EXPECT_OK(trajectory1.Append(t0, DegreesOfFreedom<ICRS>(q0, p0)));

  auto instance = ephemeris->NewInstance(
      {&trajectory1},
      Ephemeris<ICRS>::NoIntrinsicAccelerations,
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              Quinlan1999Order8A,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          10 * Second));
  EXPECT_OK(ephemeris->FlowWithFixedStep(t1, *instance));
  EXPECT_EQ(1162, trajectory1.size());

  serialization::DiscreteTrajectory message;
  trajectory1.WriteToMessage(&message, /*tracked=*/{}, /*exact=*/{});
  std::string uncompressed;
  message.SerializePartialToString(&uncompressed);
  EXPECT_EQ(18'698, uncompressed.size());

  std::string compressed;
  auto compressor = google::compression::NewGipfeliCompressor();
  compressor->Compress(uncompressed, &compressed);

  // We want a change detector, but the actual compressed size varies depending
  // on the exact numerical values, and therefore on the mathematical library.
  EXPECT_LE(17'151, compressed.size());
  EXPECT_GE(17'151, compressed.size());

  auto const trajectory2 =
      DiscreteTrajectory<ICRS>::ReadFromMessage(message, /*tracked=*/{});

  Length error;
  for (Instant t = t0; t < t1; t += 10 * Second) {
    error = std::max(
        error,
        (trajectory1.EvaluatePosition(t) - trajectory2.EvaluatePosition(t))
            .Norm());
  }
  EXPECT_THAT(error, IsNear(3.3_(1) * Metre));

  Logger logger(TEMP_DIR / "discrete_trajectory_compression.generated.wl",
                /*make_unique=*/false);
  logger.Set("trajectory1",
             trajectory1.begin(), trajectory1.end(), PreserveUnits);
  logger.Set("trajectory2",
             trajectory2.begin(), trajectory2.end(), PreserveUnits);
}

TEST(EphemerisTestNoFixture, Reanimator) {
  Instant const t_initial;
  Instant const t_final = t_initial + 5 * JulianYear;

  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");

  // Create an ephemeris and prolong it for a few years to make sure that we
  // have multiple checkpoints.
  auto ephemeris1 = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      /*fixed_step_parameters=*/{
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/10 * Minute});
  EXPECT_OK(ephemeris1->Prolong(t_final));
  EXPECT_EQ(t_initial, ephemeris1->t_min());
  EXPECT_LE(t_final, ephemeris1->t_max());

  // Serialize that ephemeris to a message and read it back.
  serialization::Ephemeris message;
  ephemeris1->WriteToMessage(&message);
  auto const ephemeris2 = Ephemeris<ICRS>::ReadFromMessage(
      /*desired_t_min=*/InfiniteFuture,
      message);

  // Reanimate the ephemeris that we just read.
  LOG(ERROR) << "Waiting until Herbert West is done...";
  ephemeris2->AwaitReanimation(t_initial);
  LOG(ERROR) << "Herbert West is finally done.";
  EXPECT_OK(ephemeris2->Prolong(t_final));

  // Check that the two ephemerides have the exact same trajectories.
  EXPECT_EQ(ephemeris1->t_min(), ephemeris2->t_min());
  EXPECT_EQ(ephemeris1->t_max(), ephemeris2->t_max());
  EXPECT_EQ(ephemeris1->bodies().size(), ephemeris2->bodies().size());
  for (int i = 0; i < ephemeris1->bodies().size(); ++i) {
    EXPECT_EQ(ephemeris1->bodies()[i]->name(), ephemeris2->bodies()[i]->name());
  }
  for (int i = 0; i < ephemeris1->bodies().size(); ++i) {
    auto trajectory1 = ephemeris1->trajectory(ephemeris1->bodies()[i]);
    auto trajectory2 = ephemeris2->trajectory(ephemeris2->bodies()[i]);
    for (Instant t = t_initial;
         t <= t_final;
         t += (t_final - t_initial) / 100) {
      EXPECT_EQ(trajectory1->EvaluateDegreesOfFreedom(t),
                trajectory2->EvaluateDegreesOfFreedom(t));
    }
  }
}
#endif

INSTANTIATE_TEST_SUITE_P(
    AllEphemerisTests,
    EphemerisTest,
    ::testing::Values(&SymplecticRungeKuttaNyströmIntegrator<
                          McLachlanAtela1992Order5Optimal,
                          Ephemeris<ICRS>::NewtonianMotionEquation>(),
                      &SymmetricLinearMultistepIntegrator<
                          Quinlan1999Order8A,
                          Ephemeris<ICRS>::NewtonianMotionEquation>()));

}  // namespace physics
}  // namespace principia
