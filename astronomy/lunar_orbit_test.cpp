
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
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {

using astronomy::ICRS;
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
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::IsNear;
using testing_utilities::PearsonProductMomentCorrelationCoefficient;
using testing_utilities::RelativeError;
using testing_utilities::Slope;

namespace astronomy {

class LunarOrbitTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    ephemeris_ = solar_system_2000_.MakeEphemeris(
        /*fitting_tolerance=*/5 * Milli(Metre),
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<ICRS>>(),
            /*step=*/10 * Minute));
  }

  static SolarSystem<ICRS> solar_system_2000_;
  static std::unique_ptr<Ephemeris<ICRS>> ephemeris_;
};

SolarSystem<ICRS> LunarOrbitTest::solar_system_2000_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2451545_000000000.proto.txt");
std::unique_ptr<Ephemeris<ICRS>> LunarOrbitTest::ephemeris_;

#if !defined(_DEBUG)

TEST_F(LunarOrbitTest, NearCircularRepeatGroundTrackOrbit) {
  auto const moon_body = dynamic_cast_not_null<OblateBody<ICRS> const*>(
      solar_system_2000_.massive_body(*ephemeris_, "Moon"));
  auto const moon_degrees_of_freedom =
      solar_system_2000_.degrees_of_freedom("Moon");

  Time const integration_step = 10 * Second;
  Time const sidereal_day = Day * 365.2425 / 366.2425;

  // Initial elements from Table 2 of Russel and Lara (2006), Repeat Ground
  // Track Lunar Orbits in the Full-Potential Plus Third-Body Problem.
  // TODO(egg): use the cartesian initial state in the MCMF frame.
  KeplerianElements<ICRS> initial_elements;
  initial_elements.semimajor_axis              = +1.861791339407e+03 *
                                                     Kilo(Metre);
  initial_elements.eccentricity                = +2.110475283361e-02;
  initial_elements.inclination                 = +9.298309294740e+01 * Degree;
  initial_elements.argument_of_periapsis       = -7.839337618501e+01 * Degree;
  initial_elements.longitude_of_ascending_node = -1.589469097527e+02 * Degree;

  MasslessBody const satellite{};
  KeplerOrbit<ICRS> initial_orbit(
      *moon_body, satellite, initial_elements, J2000);
  auto const satellite_state_vectors = initial_orbit.StateVectors(J2000);

  Time const integration_duration =
      2 * 328 * *initial_orbit.elements_at_epoch().period;

  DiscreteTrajectory<ICRS> trajectory;
  trajectory.Append(J2000, moon_degrees_of_freedom + satellite_state_vectors);
  auto const instance = ephemeris_->NewInstance(
      {&trajectory},
      Ephemeris<ICRS>::NoIntrinsicAccelerations,
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                             Position<ICRS>>(),
          integration_step));

  // Remember that because of #228 we need to loop over FlowWithFixedStep.
  ephemeris_->FlowWithFixedStep(J2000 + integration_duration, *instance);

  // Drop the units when logging to Mathematica, because it is ridiculously slow
  // at parsing them.
  base::OFStream file(SOLUTION_DIR / "mathematica" /
                      "lunar_orbit.generated.wl");
  std::vector<geometry::Vector<double, ICRS>> mma_displacements;
  std::vector<double> mma_arguments_of_periapsis;
  std::vector<double> mma_longitudes_of_ascending_nodes;

  std::vector<Angle> longitudes_of_ascending_nodes;
  std::vector<Time> times;

  for (Instant t = J2000; t <= J2000 + integration_duration;
       t += integration_duration / 100'000.0) {
    RelativeDegreesOfFreedom<ICRS> const relative_dof =
        trajectory.EvaluateDegreesOfFreedom(t) -
        ephemeris_->trajectory(moon_body)->EvaluateDegreesOfFreedom(t);
    KeplerOrbit<ICRS> actual_orbit(*moon_body, satellite, relative_dof, t);
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

    mma_displacements.push_back(relative_dof.displacement() / Metre);
    mma_arguments_of_periapsis.push_back(
        *actual_elements.argument_of_periapsis / Radian);
    mma_longitudes_of_ascending_nodes.push_back(
        actual_elements.longitude_of_ascending_node / Radian);
  }

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
