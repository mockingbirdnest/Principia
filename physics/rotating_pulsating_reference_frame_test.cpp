#include "physics/rotating_pulsating_reference_frame.hpp"

#include <memory>
#include <utility>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "physics/solar_system.hpp"
#include "testing_utilities/solar_system_factory.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace physics {

using ::testing::Eq;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_frames;
using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_generalized_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics_matchers;
using namespace principia::testing_utilities::_solar_system_factory;

class RotatingPulsatingReferenceFrameTest : public ::testing::Test {
 protected:
  using EarthMoon = Frame<struct EarthMoonTag, Arbitrary>;

  RotatingPulsatingReferenceFrameTest()
      : solar_system_(make_not_null_unique<SolarSystem<ICRS>>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt")),
        ephemeris_(solar_system_->MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            {SymmetricLinearMultistepIntegrator<
                 QuinlanTremaine1990Order12,
                 Ephemeris<ICRS>::NewtonianMotionEquation>(),
             /*step=*/10 * Minute})),
        earth_(solar_system_->massive_body(
            *ephemeris_,
            SolarSystemFactory::name(SolarSystemFactory::Earth))),
        moon_(solar_system_->massive_body(
            *ephemeris_,
            SolarSystemFactory::name(SolarSystemFactory::Moon))),
        earth_moon_(ephemeris_.get(), earth_, moon_) {}

  not_null<std::unique_ptr<SolarSystem<ICRS>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> const ephemeris_;
  not_null<MassiveBody const*> const earth_;
  not_null<MassiveBody const*> const moon_;
  _body_centred_body_direction_reference_frame::BodyCentredBodyDirectionReferenceFrame<ICRS, EarthMoon> const earth_moon_;
};

// Check that GeometricAcceleration is consistent with
// From/ToThisFrameAtTimeSimilarly.
TEST_F(RotatingPulsatingReferenceFrameTest, GeometricAcceleration) {
  Instant const t_final = J2000 + 1 * Day;
  EXPECT_OK(ephemeris_->Prolong(t_final));

  // Create trajectories in the inertial and rotating-pulsating frame,
  // from the same, starting from the same state, as well as a trajectory in the
  // inertial frame starting from a millimetrically state.
  DiscreteTrajectory<ICRS> icrs_trajectory;
  DiscreteTrajectory<ICRS> perturbed_icrs_trajectory;
  DiscreteTrajectory<EarthMoon> earth_moon_trajectory;
  EXPECT_OK(icrs_trajectory.Append(
      J2000,
      Barycentre(
          std::pair{
              ephemeris_->trajectory(earth_)->EvaluateDegreesOfFreedom(J2000),
              ephemeris_->trajectory(moon_)->EvaluateDegreesOfFreedom(J2000)},
          std::pair{1.0, 1.0})));
  EXPECT_OK(perturbed_icrs_trajectory.Append(
      icrs_trajectory.front().time,
      icrs_trajectory.front().degrees_of_freedom +
          RelativeDegreesOfFreedom<ICRS>{
              Displacement<ICRS>(
                  {1 * Milli(Metre), 1 * Milli(Metre), 1 * Milli(Metre)}),
              ICRS::unmoving}));
  EXPECT_OK(earth_moon_trajectory.Append(
      icrs_trajectory.front().time,
      earth_moon_.ToThisFrameAtTimeSimilarly(icrs_trajectory.front().time)(
          icrs_trajectory.front().degrees_of_freedom)));

  Length const length_tolerance = 1e-9 * Metre;
  Speed const speed_tolerance = 1e-9 * Metre / Second;

  // Flow the inertial trajectory.
  EXPECT_OK(ephemeris_->FlowWithAdaptiveStep(
      &icrs_trajectory,
      Ephemeris<ICRS>::NoIntrinsicAcceleration,
      t_final,
      Ephemeris<ICRS>::GeneralizedAdaptiveStepParameters(
          EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
              Fine1987RKNG34,
              Ephemeris<ICRS>::GeneralizedNewtonianMotionEquation>(),
          std::numeric_limits<std::int64_t>::max(),
          length_tolerance,
          speed_tolerance),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps));
  EXPECT_OK(ephemeris_->FlowWithAdaptiveStep(
      &perturbed_icrs_trajectory,
      Ephemeris<ICRS>::NoIntrinsicAcceleration,
      t_final,
      Ephemeris<ICRS>::GeneralizedAdaptiveStepParameters(
          EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
              Fine1987RKNG34,
              Ephemeris<ICRS>::GeneralizedNewtonianMotionEquation>(),
          std::numeric_limits<std::int64_t>::max(),
          length_tolerance,
          speed_tolerance),
      Ephemeris<ICRS>::unlimited_max_ephemeris_steps));

  // Flow the rotating-pulsating trajectory, using the geometric acceleration.
  using ODE =
      ExplicitSecondOrderOrdinaryDifferentialEquation<Position<EarthMoon>>;
  InitialValueProblem<ODE> problem{
      .equation{
          .compute_acceleration =
              [this](
                  Instant const& t,
                  std::vector<Position<EarthMoon>> const& positions,
                  std::vector<Velocity<EarthMoon>> const& velocities,
                  std::vector<Vector<Acceleration, EarthMoon>>& accelerations) {
                accelerations.front() = earth_moon_.GeometricAcceleration(
                    t, {positions.front(), velocities.front()});
                return absl::OkStatus();
              }},
      .initial_state{
          earth_moon_trajectory.front().time,
          {earth_moon_trajectory.front().degrees_of_freedom.position()},
          {earth_moon_trajectory.front().degrees_of_freedom.velocity()}}};
  auto const instance =
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Fine1987RKNG34,
                                                             ODE>()
          .NewInstance(
              problem,
              /*append_state=*/
              [&earth_moon_trajectory](ODE::State const& state) {
                EXPECT_OK(earth_moon_trajectory.Append(
                    state.time.value,
                    {state.positions.front().value,
                     state.velocities.front().value}));
              },
              [length_tolerance, speed_tolerance](
                  Time const& current_step_size,
                  ODE::State const& /*state*/,
                  ODE::State::Error const& error) {
                return std::min(
                    length_tolerance / error.position_error.front().Norm(),
                    speed_tolerance / error.velocity_error.front().Norm());
              },
              {/*first_time_step=*/t_final - J2000,
               /*safety_factor=*/0.9,
               /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
               /*last_step_is_exact=*/true});
  EXPECT_OK(instance->Solve(t_final));

  EXPECT_THAT(icrs_trajectory.back().time, Eq(t_final));
  EXPECT_THAT(perturbed_icrs_trajectory.back().time, Eq(t_final));
  EXPECT_THAT(earth_moon_trajectory.back().time, Eq(t_final));
  EXPECT_THAT(
      perturbed_icrs_trajectory.back().degrees_of_freedom,
      Componentwise(AbsoluteErrorFrom(
                        icrs_trajectory.back().degrees_of_freedom.position(),
                        IsNear(0_(1) * Milli(Metre))),
                    AbsoluteErrorFrom(
                        icrs_trajectory.back().degrees_of_freedom.velocity(),
                        IsNear(0_(1) * Milli(Metre) / Second))));
  EXPECT_THAT(
      earth_moon_.FromThisFrameAtTimeSimilarly(t_final)(
          earth_moon_trajectory.back().degrees_of_freedom),
      Componentwise(AbsoluteErrorFrom(
                        icrs_trajectory.back().degrees_of_freedom.position(),
                        IsNear(0_(1) * Milli(Metre))),
                    AbsoluteErrorFrom(
                        icrs_trajectory.back().degrees_of_freedom.velocity(),
                        IsNear(0_(1) * Milli(Metre) / Second))));
}

}  // namespace physics
}  // namespace principia
