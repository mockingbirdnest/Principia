#include <filesystem>
#include <memory>
#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_UNICODE_PATH.
#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {
namespace astronomy {

using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_frames;
using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_oblate_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_astronomy;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics;
using namespace principia::testing_utilities::_statistics;

class МолнияOrbitTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    ephemeris_ = solar_system_2000_.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<
                QuinlanTremaine1990Order12,
                Ephemeris<ICRS>::NewtonianMotionEquation>(),
            /*step=*/10 * Minute));
  }

  static SolarSystem<ICRS> solar_system_2000_;
  static std::unique_ptr<Ephemeris<ICRS>> ephemeris_;
};

SolarSystem<ICRS> МолнияOrbitTest::solar_system_2000_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2451545_000000000.proto.txt");
std::unique_ptr<Ephemeris<ICRS>> МолнияOrbitTest::ephemeris_;

#if !defined(_DEBUG)

TEST_F(МолнияOrbitTest, DISABLED_Satellite) {
  auto const earth_body = dynamic_cast_not_null<OblateBody<ICRS> const*>(
      solar_system_2000_.massive_body(*ephemeris_, "Earth"));
  auto const earth_degrees_of_freedom =
      solar_system_2000_.degrees_of_freedom("Earth");

  Time const integration_duration = 1.0 * JulianYear;
  Time const integration_step = 10 * Second;
  Time const sidereal_day = Day * 365.2425 / 366.2425;

  // These data are from https://en.wikipedia.org/wiki/Molniya_orbit.  The
  // eccentricity is from the "External links" section.
  KeplerianElements<ICRS> initial_elements;
  initial_elements.eccentricity = 0.74105;
  initial_elements.mean_motion = 2.0 * π * Radian / (sidereal_day / 2.0);
  initial_elements.inclination = ArcSin(2.0 / Sqrt(5.0));
  initial_elements.argument_of_periapsis = -π / 2.0 * Radian;
  initial_elements.longitude_of_ascending_node = 1 * Radian;
  initial_elements.mean_anomaly = 2 * Radian;

  MasslessBody const satellite{};
  KeplerOrbit<ICRS> initial_orbit(
      *earth_body, satellite, initial_elements, J2000);
  auto const satellite_state_vectors = initial_orbit.StateVectors(J2000);

  DiscreteTrajectory<ICRS> trajectory;
  EXPECT_OK(trajectory.Append(
      J2000, earth_degrees_of_freedom + satellite_state_vectors));
  auto const instance = ephemeris_->NewInstance(
      {&trajectory},
      Ephemeris<ICRS>::NoIntrinsicAccelerations,
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              Quinlan1999Order8A,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          integration_step));

  // Remember that because of #228 we need to loop over FlowWithFixedStep.
  for (Instant t = J2000 + integration_step / 2.0;
       t <= J2000 + integration_duration;
       t += integration_step / 2.0) {
    EXPECT_OK(ephemeris_->FlowWithFixedStep(t, *instance));
  }

  Logger logger(
      SOLUTION_DIR / "mathematica" /
          PRINCIPIA_UNICODE_PATH("молния_orbit.generated.wl"),
      /*make_unique=*/false);

  std::vector<Angle> longitudes_of_ascending_nodes;
  std::vector<Time> times;

  for (Instant t = J2000; t <= J2000 + integration_duration;
       t += integration_duration / 100000.0) {
    RelativeDegreesOfFreedom<ICRS> const relative_dof =
        trajectory.EvaluateDegreesOfFreedom(t) -
        ephemeris_->trajectory(earth_body)->EvaluateDegreesOfFreedom(t);
    KeplerOrbit<ICRS> actual_orbit(*earth_body, satellite, relative_dof, t);
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

    logger.Append("ppaDisplacements",
                  relative_dof.displacement(),
                  ExpressIn(Metre));
    logger.Append("ppaArguments",
                  *actual_elements.argument_of_periapsis,
                  ExpressIn(Radian));
    logger.Append("ppaLongitudes",
                  actual_elements.longitude_of_ascending_node,
                  ExpressIn(Radian));
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
  EXPECT_THAT(RelativeError(ΔΩ_per_period / (sidereal_day / 2.0),
                            actual_precession_speed),
              IsNear(0.076_(1)));
}

#endif

}  // namespace astronomy
}  // namespace principia
