
#include <numeric>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using base::not_null;
using geometry::AngleBetween;
using integrators::McLachlanAtela1992Order5Optimal;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::RelativeDegreesOfFreedom;
using physics::SolarSystem;
using quantities::astronomy::JulianYear;
using quantities::si::AstronomicalUnit;
using quantities::si::Day;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;
using testing_utilities::AbsoluteError;

namespace astronomy {

class MercuryPerihelionTest : public testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    solar_system_1950_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2433282_500000000.proto.txt");
    ephemeris_ = solar_system_1950_.MakeEphemeris(
        /*fitting_tolerance=*/5 * Milli(Metre),
        Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
            McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
            /*step=*/45 * Minute));
  }

  MercuryPerihelionTest()
      : sun_(solar_system_1950_.massive_body(*ephemeris_, "Sun")),
        mercury_(solar_system_1950_.massive_body(*ephemeris_, "Mercury")) {}

  static SolarSystem<ICRFJ2000Equator> solar_system_1950_;
  static std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
  not_null<MassiveBody const*> sun_;
  not_null<MassiveBody const*> mercury_;
};

SolarSystem<ICRFJ2000Equator> MercuryPerihelionTest::solar_system_1950_;
std::unique_ptr<Ephemeris<ICRFJ2000Equator>> MercuryPerihelionTest::ephemeris_;

TEST_F(MercuryPerihelionTest, PrintPerihelion) {
  // From Horizons by dichotomy, around 1950-JAN-11 03:12:30.0000 (TDB).
  Instant const first_perihelion_low = JulianDate(2433292.633686343);
  Instant const first_perihelion_high = JulianDate(2433292.633692130);
  Instant const first_perihelion_mid = Barycentre<Instant, double>(
      {first_perihelion_low, first_perihelion_high}, {1, 1});

  // From Horizons by dichotomy, around 1959-Nov-26 21:00:28.0000 (TDB).
  Instant const last_perihelion_low = JulianDate(2436899.375324074);
  Instant const last_perihelion_high = JulianDate(2436899.375329861 );
  Instant const last_perihelion_mid = Barycentre<Instant, double>(
      {last_perihelion_low, last_perihelion_high}, {1, 1});

  DiscreteTrajectory<ICRFJ2000Equator> sun_apoapsides;
  DiscreteTrajectory<ICRFJ2000Equator> sun_periapsides;
  DiscreteTrajectory<ICRFJ2000Equator> mercury_apoapsides;
  DiscreteTrajectory<ICRFJ2000Equator> mercury_periapsides;
  ephemeris_->Prolong(solar_system_1950_.epoch() + 10 * JulianYear);
  ephemeris_->ComputeApsides(sun_,
                             mercury_,
                             sun_apoapsides,
                             sun_periapsides,
                             mercury_apoapsides,
                             mercury_periapsides);

  EXPECT_LT(AbsoluteError(sun_periapsides.Begin().time(), first_perihelion_mid),
            0.26 * Second);
  EXPECT_LT(AbsoluteError(sun_periapsides.last().time(), last_perihelion_mid),
            99.3 * Second);

  std::experimental::optional<Instant> previous_time;
  std::experimental::optional<KeplerianElements<ICRFJ2000Equator>>
      previous_keplerian_elements;
  std::vector<AngularFrequency> precessions;
  for (auto sun_it = sun_periapsides.Begin(),
            mercury_it = mercury_periapsides.Begin();
       sun_it != sun_periapsides.End() &&
       mercury_it != mercury_periapsides.End();
       ++sun_it, ++mercury_it) {
    Instant const time = sun_it.time();
    RelativeDegreesOfFreedom<ICRFJ2000Equator> const
        relative_degrees_of_freedom =
            sun_it.degrees_of_freedom() - mercury_it.degrees_of_freedom();
    KeplerOrbit<ICRFJ2000Equator> orbit(
        *sun_, *mercury_, relative_degrees_of_freedom, time);
    KeplerianElements<ICRFJ2000Equator> const keplerian_elements =
        orbit.elements_at_epoch();
    Angle const anomaly =
        keplerian_elements.mean_anomaly < 1 * Radian
            ? keplerian_elements.mean_anomaly
            : 2 * π * Radian - keplerian_elements.mean_anomaly;
    EXPECT_LT(anomaly, 1e-12 * Radian);
    if (previous_time) {
      AngularFrequency const precession =
          (keplerian_elements.argument_of_periapsis -
           previous_keplerian_elements->argument_of_periapsis) /
          (time - *previous_time);
      precessions.push_back(precession);
      LOG(ERROR)<<precession * 100 * JulianYear / (1 * ArcSecond);
    }
    previous_time = time;
    previous_keplerian_elements = keplerian_elements;
  }
  AngularFrequency average;
  for (auto const precession : precessions) {
    average += precession;
  }
  average /= precessions.size();
  LOG(ERROR) << "Average: " << average * 100 * JulianYear / (1 * ArcSecond);
}

}  // namespace astronomy
}  // namespace principia
