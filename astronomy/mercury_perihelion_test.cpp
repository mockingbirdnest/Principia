
#include <fstream>
#include <numeric>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "mathematica/mathematica.hpp"
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
using physics::ContinuousTrajectory;
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
using testing_utilities::RelativeError;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;

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
        mercury_(solar_system_1950_.massive_body(*ephemeris_, "Mercury")),
        t_1950_(JulianDate(2433282.500000000)),
        t_1960_(JulianDate(2436934.500000000)),
        t_2050_(JulianDate(2469807.500000000)) {
    keplerian_elements_1950_.eccentricity = 2.056187274905493E-01;
    keplerian_elements_1950_.semimajor_axis =
        5.790897350196702E+07 * Kilo(Metre);
    keplerian_elements_1950_.mean_motion =
        4.736523381721572E-05 * Degree / Second;
    keplerian_elements_1950_.inclination = 2.854970888858858E+01 * Degree;
    keplerian_elements_1950_.longitude_of_ascending_node =
        1.100421919049157E+01 * Degree;
    keplerian_elements_1950_.argument_of_periapsis =
        6.747518782664667E+01 * Degree;
    keplerian_elements_1950_.mean_anomaly = 3.185292722214373E+02 * Degree;

    keplerian_elements_1960_.eccentricity = 2.056305163902087E-01;
    keplerian_elements_1960_.semimajor_axis =
        5.790908221204869E+07 * Kilo(Metre);
    keplerian_elements_1960_.mean_motion =
        4.736510044238512E-05 * Degree / Second;
    keplerian_elements_1960_.inclination = 2.855025264340049E+01 * Degree;
    keplerian_elements_1960_.longitude_of_ascending_node =
        1.100120089390144E+01 * Degree;
    keplerian_elements_1960_.argument_of_periapsis =
        6.748880915164143E+01 * Degree;
    keplerian_elements_1960_.mean_anomaly = 1.437427180144607E+02 * Degree;

    keplerian_elements_2050_.eccentricity = 2.056443534543679E-01;
    keplerian_elements_2050_.semimajor_axis =
        5.790909484661692E+07 * Kilo(Metre);
    keplerian_elements_2050_.mean_motion =
        4.736508494125626E-05 * Degree / Second;
    keplerian_elements_2050_.inclination = 2.855459021021251E+01 * Degree;
    keplerian_elements_2050_.longitude_of_ascending_node =
        1.097153493770216E+01 * Degree;
    keplerian_elements_2050_.argument_of_periapsis =
        6.765993779473686E+01 * Degree;
    keplerian_elements_2050_.mean_anomaly = 3.105229688140852E+01 * Degree;
  }

  static SolarSystem<ICRFJ2000Equator> solar_system_1950_;
  static std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;

  not_null<MassiveBody const*> sun_;
  not_null<MassiveBody const*> mercury_;
  Instant t_1950_;
  Instant t_1960_;
  Instant t_2050_;
  KeplerianElements<ICRFJ2000Equator> keplerian_elements_1950_;
  KeplerianElements<ICRFJ2000Equator> keplerian_elements_1960_;
  KeplerianElements<ICRFJ2000Equator> keplerian_elements_2050_;
};

SolarSystem<ICRFJ2000Equator> MercuryPerihelionTest::solar_system_1950_;
std::unique_ptr<Ephemeris<ICRFJ2000Equator>> MercuryPerihelionTest::ephemeris_;

TEST_F(MercuryPerihelionTest, Year1950) {
  ephemeris_->Prolong(t_1950_);

  auto const& sun_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Sun");
  auto const& mercury_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Mercury");

  RelativeDegreesOfFreedom<ICRFJ2000Equator> const relative_degrees_of_freedom =
      mercury_trajectory.EvaluateDegreesOfFreedom(t_1950_, /*hint=*/nullptr) -
      sun_trajectory.EvaluateDegreesOfFreedom(t_1950_, /*hint=*/nullptr);
  KeplerOrbit<ICRFJ2000Equator> orbit(
      *sun_, *mercury_, relative_degrees_of_freedom, t_1950_);
  KeplerianElements<ICRFJ2000Equator> const keplerian_elements =
      orbit.elements_at_epoch();

  EXPECT_LT(RelativeError(keplerian_elements.eccentricity,
                          keplerian_elements_1950_.eccentricity), 3.8e-14);
  EXPECT_LT(RelativeError(*keplerian_elements.semimajor_axis,
                          *keplerian_elements_1950_.semimajor_axis), 1.1e-14);
  EXPECT_LT(RelativeError(*keplerian_elements.mean_motion,
                          *keplerian_elements_1950_.mean_motion), 1.5e-14);
  EXPECT_LT(RelativeError(keplerian_elements.inclination,
                          keplerian_elements_1950_.inclination), 5.2e-15);
  EXPECT_LT(RelativeError(keplerian_elements.longitude_of_ascending_node,
                          keplerian_elements_1950_.
                              longitude_of_ascending_node), 4.2e-15);
  EXPECT_LT(RelativeError(keplerian_elements_1950_.argument_of_periapsis,
                          keplerian_elements.argument_of_periapsis), 8.5e-14);
  EXPECT_LT(RelativeError(keplerian_elements.mean_anomaly,
                          keplerian_elements_1950_.mean_anomaly), 1.3e-14);
}

TEST_F(MercuryPerihelionTest, Year1960) {
  ephemeris_->Prolong(t_1960_);

  auto const& sun_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Sun");
  auto const& mercury_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Mercury");

  RelativeDegreesOfFreedom<ICRFJ2000Equator> const relative_degrees_of_freedom =
      mercury_trajectory.EvaluateDegreesOfFreedom(t_1960_, /*hint=*/nullptr) -
      sun_trajectory.EvaluateDegreesOfFreedom(t_1960_, /*hint=*/nullptr);
  KeplerOrbit<ICRFJ2000Equator> orbit(
      *sun_, *mercury_, relative_degrees_of_freedom, t_1960_);
  KeplerianElements<ICRFJ2000Equator> const keplerian_elements =
      orbit.elements_at_epoch();

  EXPECT_LT(RelativeError(keplerian_elements.eccentricity,
                          keplerian_elements_1960_.eccentricity), 5.3e-7);
  EXPECT_LT(RelativeError(*keplerian_elements.semimajor_axis,
                          *keplerian_elements_1960_.semimajor_axis), 1.1e-7);
  EXPECT_LT(RelativeError(*keplerian_elements.mean_motion,
                          *keplerian_elements_1960_.mean_motion), 1.7e-7);
  EXPECT_LT(RelativeError(keplerian_elements.inclination,
                          keplerian_elements_1960_.inclination), 1.1e-10);
  EXPECT_LT(RelativeError(keplerian_elements.longitude_of_ascending_node,
                          keplerian_elements_1960_.longitude_of_ascending_node),
            1.3e-9);
  EXPECT_THAT(keplerian_elements_1960_.argument_of_periapsis -
                  keplerian_elements.argument_of_periapsis,
              AllOf(Gt(4.202 * ArcSecond), Lt(4.203 * ArcSecond)));
  EXPECT_THAT(keplerian_elements.mean_anomaly -
                  keplerian_elements_1960_.mean_anomaly,
              AllOf(Gt(17.0 * ArcSecond), Lt(17.1 * ArcSecond)));
}

TEST_F(MercuryPerihelionTest, KeplerianElements) {
  LOG(ERROR) << JulianDate(2469807.500000000);
  ephemeris_->Prolong(solar_system_1950_.epoch() + 100 * JulianYear);

  auto const& sun_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Sun");
  auto const& mercury_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Mercury");
  typename ContinuousTrajectory<ICRFJ2000Equator>::Hint sun_hint;
  typename ContinuousTrajectory<ICRFJ2000Equator>::Hint mercury_hint;

  std::vector<Angle> arguments_of_periapsis;
  std::experimental::optional<Instant> previous_time;
  std::experimental::optional<KeplerianElements<ICRFJ2000Equator>>
      previous_keplerian_elements;
  for (Instant time = solar_system_1950_.epoch();
       time <= solar_system_1950_.epoch() + 100 * JulianYear;
       time += 1 * Day) {
    RelativeDegreesOfFreedom<ICRFJ2000Equator> const
        relative_degrees_of_freedom =
            mercury_trajectory.EvaluateDegreesOfFreedom(time, &mercury_hint) -
            sun_trajectory.EvaluateDegreesOfFreedom(time, &sun_hint);
    KeplerOrbit<ICRFJ2000Equator> orbit(
        *sun_, *mercury_, relative_degrees_of_freedom, time);
    KeplerianElements<ICRFJ2000Equator> const keplerian_elements =
        orbit.elements_at_epoch();
    arguments_of_periapsis.push_back(keplerian_elements.argument_of_periapsis);
    if (!previous_keplerian_elements) {
      LOG(ERROR)<<time<<":\n"<<keplerian_elements;
    }
    previous_time = time;
    previous_keplerian_elements = keplerian_elements;
  }
  LOG(ERROR) << *previous_time << ":\n" << *previous_keplerian_elements;

  std::ofstream file;
  file.open("mercury_perihelion.generated.wl");
  file << mathematica::Assign("argumentsOfPeriapsis", arguments_of_periapsis);
  file.close();
}

TEST_F(MercuryPerihelionTest, PrintPerihelion) {
  // From Horizons by dichotomy, around 1950-01-11T03:12:30.0000Z (TDB).
  Instant const first_perihelion_time_low = JulianDate(2433292.633686343);
  Instant const first_perihelion_time_high = JulianDate(2433292.633692130);
  Instant const first_perihelion_time_mid = Barycentre<Instant, double>(
      {first_perihelion_time_low, first_perihelion_time_high}, {1, 1});

  // From Horizons by dichotomy, around 1959-11-26T21:00:28.0000Z (TDB).
  Instant const last_perihelion_time_low = JulianDate(2436899.375324074);
  Instant const last_perihelion_time_high = JulianDate(2436899.375329861 );
  Instant const last_perihelion_time_mid = Barycentre<Instant, double>(
      {last_perihelion_time_low, last_perihelion_time_high}, {1, 1});

  Position<ICRFJ2000Equator> const last_perihelion_position_low =
      ICRFJ2000Equator::origin +
      Displacement<ICRFJ2000Equator>(
          {7.084639885659656e-02 * AstronomicalUnit,
           2.748967000176096e-01 * AstronomicalUnit,
           1.388562317323265e-01 * AstronomicalUnit});
  Position<ICRFJ2000Equator> const last_perihelion_position_high =
      ICRFJ2000Equator::origin +
      Displacement<ICRFJ2000Equator>(
          {7.084620740779489e-02 * AstronomicalUnit,
           2.748967303484739e-01 * AstronomicalUnit,
           1.388562678052315e-01 * AstronomicalUnit});
  Position<ICRFJ2000Equator> const last_perihelion_position_mid =
      Barycentre<Position<ICRFJ2000Equator>, double>(
          {last_perihelion_position_low, last_perihelion_position_high},
          {1, 1});

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

  EXPECT_LT(AbsoluteError(sun_periapsides.Begin().time(),
                          first_perihelion_time_mid),
            0.26 * Second);
  EXPECT_LT(AbsoluteError(sun_periapsides.last().time(),
                          last_perihelion_time_mid),
            99.3 * Second);
  EXPECT_LT(
      AbsoluteError(mercury_periapsides.last().degrees_of_freedom().position(),
                    last_perihelion_position_mid),
      946 * Kilo(Metre));

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
