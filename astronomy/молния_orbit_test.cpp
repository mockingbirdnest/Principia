
#include <experimental/filesystem>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using astronomy::J2000;
using geometry::Position;
using integrators::Quinlan1999Order8A;
using integrators::QuinlanTremaine1990Order12;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::KeplerOrbit;
using physics::KeplerianElements;
using physics::MasslessBody;
using physics::SolarSystem;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::ArcSin;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;

namespace astronomy {

class МолнияOrbitTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    solar_system_2000_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2451545_000000000.proto.txt");
    ephemeris_ = solar_system_2000_.MakeEphemeris(
        /*fitting_tolerance=*/5 * Milli(Metre),
        Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
            QuinlanTremaine1990Order12<Position<ICRFJ2000Equator>>(),
            /*step=*/10 * Minute));
  }

  static SolarSystem<ICRFJ2000Equator> solar_system_2000_;
  static std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
};

SolarSystem<ICRFJ2000Equator> МолнияOrbitTest::solar_system_2000_;
std::unique_ptr<Ephemeris<ICRFJ2000Equator>> МолнияOrbitTest::ephemeris_;

TEST_F(МолнияOrbitTest, Satellite) {
  auto const earth_body = solar_system_2000_.massive_body(*ephemeris_, "Earth");
  auto const earth_degrees_of_freedom =
      solar_system_2000_.initial_state("Earth");

  Time const sidereal_day = Day * 365.2425 / 366.2425;

  // These data are from https://en.wikipedia.org/wiki/Molniya_orbit.  The
  // eccentricity is from the "External links" section.
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity = 0.74105;
  elements.mean_motion = 2.0 * π * Radian / (sidereal_day / 2.0);
  elements.inclination = ArcSin(2.0 / Sqrt(5.0));
  elements.longitude_of_ascending_node = -π / 2.0 * Radian;
  elements.longitude_of_ascending_node = 1.0 * Radian;
  elements.mean_anomaly = 2.0 * Radian;

  MasslessBody const satellite;
  KeplerOrbit<ICRFJ2000Equator> orbit(*earth_body, satellite, elements, J2000);
  auto const satellite_state_vectors = orbit.StateVectors(J2000);

  DiscreteTrajectory<ICRFJ2000Equator> trajectory;
  trajectory.Append(J2000, earth_degrees_of_freedom + satellite_state_vectors);
  auto const instance = ephemeris_->NewInstance(
      {&trajectory},
      Ephemeris<ICRFJ2000Equator>::NoIntrinsicAccelerations,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
          Quinlan1999Order8A<Position<ICRFJ2000Equator>>(),
          10 * Second));

  ephemeris_->FlowWithFixedStep(J2000 + 1.0 * JulianYear, *instance);
}

}  // namespace astronomy
}  // namespace principia
