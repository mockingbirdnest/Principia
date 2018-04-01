
#include <filesystem>
#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "base/file.hpp"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using astronomy::J2000;
using base::dynamic_cast_not_null;
using base::OFStream;
using geometry::Instant;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::KeplerOrbit;
using physics::KeplerianElements;
using physics::MasslessBody;
using physics::OblateBody;
using physics::RelativeDegreesOfFreedom;
using physics::SolarSystem;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::ArcSin;
using quantities::Cos;
using quantities::Length;
using quantities::Pow;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::PearsonProductMomentCorrelationCoefficient;
using testing_utilities::RelativeError;
using testing_utilities::Slope;

namespace astronomy {

class МолнияOrbitTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    ephemeris_ = solar_system_2000_.MakeEphemeris(
        /*fitting_tolerance=*/5 * Milli(Metre),
        Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<ICRFJ2000Equator>>(),
            /*step=*/10 * Minute));
  }

  static SolarSystem<ICRFJ2000Equator> solar_system_2000_;
  static std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
};

SolarSystem<ICRFJ2000Equator> МолнияOrbitTest::solar_system_2000_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2451545_000000000.proto.txt");
std::unique_ptr<Ephemeris<ICRFJ2000Equator>> МолнияOrbitTest::ephemeris_;

#if !defined(_DEBUG)

TEST_F(МолнияOrbitTest, Satellite) {
  auto const earth_body =
      dynamic_cast_not_null<OblateBody<ICRFJ2000Equator> const*>(
          solar_system_2000_.massive_body(*ephemeris_, "Earth"));
  auto const earth_degrees_of_freedom =
      solar_system_2000_.degrees_of_freedom("Earth");

  Time const integration_duration = 1.0 * JulianYear;
  Time const integration_step = 10 * Second;
  Time const sidereal_day = Day * 365.2425 / 366.2425;

  // These data are from https://en.wikipedia.org/wiki/Molniya_orbit.  The
  // eccentricity is from the "External links" section.
  KeplerianElements<ICRFJ2000Equator> initial_elements;
  initial_elements.eccentricity = 0.74105;
  initial_elements.mean_motion = 2.0 * π * Radian / (sidereal_day / 2.0);
  initial_elements.inclination = ArcSin(2.0 / Sqrt(5.0));
  initial_elements.argument_of_periapsis = -π / 2.0 * Radian;
  initial_elements.longitude_of_ascending_node = 1 * Radian;
  initial_elements.mean_anomaly = 2 * Radian;

  MasslessBody const satellite{};
  KeplerOrbit<ICRFJ2000Equator> initial_orbit(
      *earth_body, satellite, initial_elements, J2000);
  auto const satellite_state_vectors = initial_orbit.StateVectors(J2000);

  DiscreteTrajectory<ICRFJ2000Equator> trajectory;
  trajectory.Append(J2000, earth_degrees_of_freedom + satellite_state_vectors);
  auto const instance = ephemeris_->NewInstance(
      {&trajectory},
      Ephemeris<ICRFJ2000Equator>::NoIntrinsicAccelerations,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                             Position<ICRFJ2000Equator>>(),
          integration_step));

  // Remember that because of #228 we need to loop over FlowWithFixedStep.
  for (Instant t = J2000 + integration_step / 2.0;
       t <= J2000 + integration_duration;
       t += integration_step / 2.0) {
    ephemeris_->FlowWithFixedStep(t, *instance);
  }

  // Drop the units when logging to Mathematica, because it is ridiculously slow
  // at parsing them.
  base::OFStream file(SOLUTION_DIR / "mathematica" /
                      UNICODE_PATH("молния_orbit.generated.wl"));
  std::vector<geometry::Vector<double, ICRFJ2000Equator>> mma_displacements;
  std::vector<double> mma_arguments_of_periapsis;
  std::vector<double> mma_longitudes_of_ascending_nodes;

  std::vector<Angle> longitudes_of_ascending_nodes;
  std::vector<Time> times;

  for (Instant t = J2000; t <= J2000 + integration_duration;
       t += integration_duration / 100000.0) {
    RelativeDegreesOfFreedom<ICRFJ2000Equator> const relative_dof =
        trajectory.EvaluateDegreesOfFreedom(t) -
        ephemeris_->trajectory(earth_body)->EvaluateDegreesOfFreedom(t);
    KeplerOrbit<ICRFJ2000Equator> actual_orbit(*earth_body,
                                               satellite,
                                               relative_dof,
                                               t);
    auto actual_elements = actual_orbit.elements_at_epoch();

    if (actual_elements.longitude_of_ascending_node >
        initial_elements.longitude_of_ascending_node + π * Radian) {
      actual_elements.longitude_of_ascending_node -= 2.0 * π * Radian;
    }
    if (actual_elements.longitude_of_ascending_node <
        initial_elements.longitude_of_ascending_node - π * Radian) {
      actual_elements.longitude_of_ascending_node += 2.0 * π * Radian;
    }
    longitudes_of_ascending_nodes.push_back(
        actual_elements.longitude_of_ascending_node -
        initial_elements.longitude_of_ascending_node);
    times.push_back(t - J2000);

    // Check that the argument of the perigee remains roughly constant (modulo
    // the influence of the Moon).
    EXPECT_LT(RelativeError(
                  2.0 * π * Radian + *initial_elements.argument_of_periapsis,
                  *actual_elements.argument_of_periapsis),
              0.0026);

    mma_displacements.push_back(relative_dof.displacement() / Metre);
    mma_arguments_of_periapsis.push_back(
        *actual_elements.argument_of_periapsis / Radian);
    mma_longitudes_of_ascending_nodes.push_back(
        actual_elements.longitude_of_ascending_node / Radian);
  }

  // Check that we have a regular precession of the longitude.
  double const correlation_coefficients =
      PearsonProductMomentCorrelationCoefficient(times,
                                                 longitudes_of_ascending_nodes);
  EXPECT_GT(correlation_coefficients, -0.99999);
  EXPECT_LT(correlation_coefficients, -0.99998);

  // Check that the longitude precesses at the right speed, mostly.
  AngularFrequency const actual_precession_speed =
      Slope(times, longitudes_of_ascending_nodes);
  Length const& semilatus_rectum =
      *initial_orbit.elements_at_epoch().semilatus_rectum;
  Angle const ΔΩ_per_period = -2.0 * π * Radian * earth_body->j2_over_μ() /
                              (semilatus_rectum * semilatus_rectum) *
                              (3.0 / 2.0) * Cos(initial_elements.inclination);
  EXPECT_LT(RelativeError(ΔΩ_per_period / (sidereal_day / 2.0),
                          actual_precession_speed),
            0.076);

  file << mathematica::Assign("ppaDisplacements",
                              mma_displacements);
  file << mathematica::Assign("ppaArguments",
                              mma_arguments_of_periapsis);
  file << mathematica::Assign("ppaLongitudes",
                              mma_longitudes_of_ascending_nodes);
}

#endif

}  // namespace astronomy
}  // namespace principia
