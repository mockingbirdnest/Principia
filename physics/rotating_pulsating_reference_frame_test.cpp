#include "physics/rotating_pulsating_reference_frame.hpp"

#include <memory>
#include <utility>

#include "astronomy/orbital_elements.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "physics/body_centred_body_direction_reference_frame.hpp"
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
using namespace principia::integrators::_explicit_runge_kutta_integrator;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::physics::_body_centred_body_direction_reference_frame;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
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
  // REMOVE BEFORE FLIGHT
  BodyCentredBodyDirectionReferenceFrame<ICRS, EarthMoon> const earth_moon_;
};

#if !defined(_DEBUG)

// Check that GeometricAcceleration is consistent with
// From/ToThisFrameAtTimeSimilarly.
TEST_F(RotatingPulsatingReferenceFrameTest, GeometricAcceleration) {
  Instant const t_final = J2000 + 1 * Day;
  EXPECT_OK(ephemeris_->Prolong(t_final));

  // Create trajectories in the inertial and rotating-pulsating frame, starting
  // from the same state.
  DiscreteTrajectory<ICRS> icrs_trajectory;
  DiscreteTrajectory<EarthMoon> earth_moon_trajectory;
  EXPECT_OK(icrs_trajectory.Append(
      J2000,
      ephemeris_->trajectory(earth_)->EvaluateDegreesOfFreedom(J2000) +
          KeplerOrbit<ICRS>(
              *earth_,
              MasslessBody{},
              KeplerianElements<ICRS>{.eccentricity = 0.3,
                                      .period = 1 * Day,
                                      .inclination = 0 * Degree,
                                      .longitude_of_ascending_node = 0 * Degree,
                                      .argument_of_periapsis = 0 * Degree,
                                      .mean_anomaly = 0 * Degree},
              J2000)
              .StateVectors(J2000)));
  EXPECT_OK(earth_moon_trajectory.Append(
      icrs_trajectory.front().time,
      earth_moon_.ToThisFrameAtTimeSimilarly(icrs_trajectory.front().time)(
          icrs_trajectory.front().degrees_of_freedom)));

  auto const icrs_instance = ephemeris_->NewInstance(
      {&icrs_trajectory},
      Ephemeris<ICRS>::NoIntrinsicAccelerations,
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              Quinlan1999Order8A,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/1.0 * Second));
  // Flow the inertial trajectory.
  EXPECT_OK(ephemeris_->FlowWithFixedStep(t_final, *icrs_instance));

  // Flow the rotating-pulsating trajectory, using the geometric acceleration.
  using ODE =
      ExplicitFirstOrderOrdinaryDifferentialEquation<Instant,
                                                     Position<EarthMoon>,
                                                     Velocity<EarthMoon>>;
  InitialValueProblem<ODE> problem{
      .equation{.compute_derivative =
                    [this](Instant const& t,
                           ODE::DependentVariables const& y,
                           ODE::DependentVariableDerivatives& yʹ) {
                      auto const& [q, v] = y;
                      auto& [qʹ, vʹ] = yʹ;
                      qʹ = v;
                      vʹ = earth_moon_.GeometricAcceleration(t, {q, v});
                      return absl::OkStatus();
                    }},
      .initial_state{
          earth_moon_trajectory.front().time,
          {earth_moon_trajectory.front().degrees_of_freedom.position(),
           earth_moon_trajectory.front().degrees_of_freedom.velocity()}}};
  auto const earth_moon_instance =
      ExplicitRungeKuttaIntegrator<Kutta1901Vσ1, ODE>()
          .NewInstance(
              problem,
              /*append_state=*/
              [&earth_moon_trajectory](ODE::State const& state) {
                auto const& [q, v] = state.y;
                EXPECT_OK(earth_moon_trajectory.Append(state.s.value,
                                                       {q.value, v.value}));
              },
              /*step=*/1.0 * Second);
  EXPECT_OK(earth_moon_instance->Solve(t_final));

  EXPECT_THAT(icrs_trajectory.back().time, Eq(t_final));
  EXPECT_THAT(earth_moon_trajectory.back().time, Eq(t_final));
  // TODO(egg): These errors make no sense.
  EXPECT_THAT(
      earth_moon_.FromThisFrameAtTimeSimilarly(t_final)(
          earth_moon_trajectory.back().degrees_of_freedom),
      Componentwise(AbsoluteErrorFrom(
                        icrs_trajectory.back().degrees_of_freedom.position(),
                        IsNear(878.29_(1) * Metre)),
                    AbsoluteErrorFrom(
                        icrs_trajectory.back().degrees_of_freedom.velocity(),
                        IsNear(96.08_(1) * Milli(Metre) / Second))));
}

#endif

}  // namespace physics
}  // namespace principia
