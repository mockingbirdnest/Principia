
#include <fstream>
#include <numeric>

#include "astronomy/epoch.hpp"
#include "astronomy/time_scales.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
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
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using base::not_null;
using geometry::AngleBetween;
using geometry::Instant;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::ContinuousTrajectory;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::RelativeDegreesOfFreedom;
using physics::SolarSystem;
using quantities::astronomy::JulianYear;
using quantities::astronomy::AstronomicalUnit;
using quantities::si::ArcSecond;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Minute;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;

namespace astronomy {

class MercuryPerihelionTest : public testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    ephemeris_ = solar_system_1950_.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<ICRS>>(),
            /*step=*/10 * Minute));
  }

  MercuryPerihelionTest()
      : sun_(solar_system_1950_.massive_body(*ephemeris_, "Sun")),
        mercury_(solar_system_1950_.massive_body(*ephemeris_, "Mercury")),
        t_1950_("1950-01-01T00:00:00"_TT),
        t_1960_("1960-01-01T00:00:00"_TT),
        t_2050_("2050-01-01T00:00:00"_TT) {
    keplerian_elements_1950_.eccentricity = 2.056187274905493e-01;
    keplerian_elements_1950_.semimajor_axis =
        5.790897350196702e+07 * Kilo(Metre);
    keplerian_elements_1950_.mean_motion =
        4.736523381721572e-05 * Degree / Second;
    keplerian_elements_1950_.inclination = 2.854970888858858e+01 * Degree;
    keplerian_elements_1950_.longitude_of_ascending_node =
        1.100421919049157e+01 * Degree;
    keplerian_elements_1950_.argument_of_periapsis =
        6.747518782664667e+01 * Degree;
    keplerian_elements_1950_.mean_anomaly = 3.185292722214373e+02 * Degree;

    keplerian_elements_1960_.eccentricity = 2.056305163902087e-01;
    keplerian_elements_1960_.semimajor_axis =
        5.790908221204869e+07 * Kilo(Metre);
    keplerian_elements_1960_.mean_motion =
        4.736510044238512e-05 * Degree / Second;
    keplerian_elements_1960_.inclination = 2.855025264340049e+01 * Degree;
    keplerian_elements_1960_.longitude_of_ascending_node =
        1.100120089390144e+01 * Degree;
    keplerian_elements_1960_.argument_of_periapsis =
        6.748880915164143e+01 * Degree;
    keplerian_elements_1960_.mean_anomaly = 1.437427180144607e+02 * Degree;

    keplerian_elements_2050_.eccentricity = 2.056443534543679e-01;
    keplerian_elements_2050_.semimajor_axis =
        5.790909484661692e+07 * Kilo(Metre);
    keplerian_elements_2050_.mean_motion =
        4.736508494125626e-05 * Degree / Second;
    keplerian_elements_2050_.inclination = 2.855459021021251e+01 * Degree;
    keplerian_elements_2050_.longitude_of_ascending_node =
        1.097153493770216e+01 * Degree;
    keplerian_elements_2050_.argument_of_periapsis =
        6.765993779473686e+01 * Degree;
    keplerian_elements_2050_.mean_anomaly = 3.105229688140852e+01 * Degree;
  }

  static SolarSystem<ICRS> solar_system_1950_;
  static std::unique_ptr<Ephemeris<ICRS>> ephemeris_;

  not_null<MassiveBody const*> sun_;
  not_null<MassiveBody const*> mercury_;
  Instant t_1950_;
  Instant t_1960_;
  Instant t_2050_;
  KeplerianElements<ICRS> keplerian_elements_1950_;
  KeplerianElements<ICRS> keplerian_elements_1960_;
  KeplerianElements<ICRS> keplerian_elements_2050_;
};

SolarSystem<ICRS> MercuryPerihelionTest::solar_system_1950_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2433282_500000000.proto.txt");
std::unique_ptr<Ephemeris<ICRS>> MercuryPerihelionTest::ephemeris_;

TEST_F(MercuryPerihelionTest, Year1950) {
  ephemeris_->Prolong(t_1950_);

  auto const& sun_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Sun");
  auto const& mercury_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Mercury");

  RelativeDegreesOfFreedom<ICRS> const relative_degrees_of_freedom =
      mercury_trajectory.EvaluateDegreesOfFreedom(t_1950_) -
      sun_trajectory.EvaluateDegreesOfFreedom(t_1950_);
  KeplerOrbit<ICRS> orbit(
      *sun_, *mercury_, relative_degrees_of_freedom, t_1950_);
  KeplerianElements<ICRS> const keplerian_elements = orbit.elements_at_epoch();

  EXPECT_THAT(RelativeError(*keplerian_elements.eccentricity,
                             *keplerian_elements_1950_.eccentricity),
              IsNear(3.6e-13));
  EXPECT_THAT(RelativeError(*keplerian_elements.semimajor_axis,
                            *keplerian_elements_1950_.semimajor_axis),
              IsNear(7.7e-14));
  EXPECT_THAT(RelativeError(*keplerian_elements.mean_motion,
                          *keplerian_elements_1950_.mean_motion),
              IsNear(1.2e-13));
  EXPECT_THAT(RelativeError(keplerian_elements.inclination,
                          keplerian_elements_1950_.inclination),
              IsNear(1.7e-14));
  EXPECT_THAT(RelativeError(keplerian_elements.longitude_of_ascending_node,
                           keplerian_elements_1950_.
                               longitude_of_ascending_node),
              IsNear(1.2e-14));
  EXPECT_THAT(RelativeError(*keplerian_elements_1950_.argument_of_periapsis,
                          *keplerian_elements.argument_of_periapsis),
              IsNear(6.6e-14));
  EXPECT_THAT(RelativeError(*keplerian_elements.mean_anomaly,
                          *keplerian_elements_1950_.mean_anomaly),
              IsNear(3.0e-14));
}

#if !defined(_DEBUG)
TEST_F(MercuryPerihelionTest, Year1960) {
  ephemeris_->Prolong(t_1960_);

  auto const& sun_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Sun");
  auto const& mercury_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Mercury");

  RelativeDegreesOfFreedom<ICRS> const relative_degrees_of_freedom =
      mercury_trajectory.EvaluateDegreesOfFreedom(t_1960_) -
      sun_trajectory.EvaluateDegreesOfFreedom(t_1960_);
  KeplerOrbit<ICRS> orbit(
      *sun_, *mercury_, relative_degrees_of_freedom, t_1960_);
  KeplerianElements<ICRS> const keplerian_elements = orbit.elements_at_epoch();

  EXPECT_THAT(RelativeError(*keplerian_elements.eccentricity,
                            *keplerian_elements_1960_.eccentricity),
              IsNear(5.3e-7));
  EXPECT_THAT(RelativeError(*keplerian_elements.semimajor_axis,
                            *keplerian_elements_1960_.semimajor_axis),
              IsNear(1.1e-7));
  EXPECT_THAT(RelativeError(*keplerian_elements.mean_motion,
                            *keplerian_elements_1960_.mean_motion),
              IsNear(1.6e-7));
  EXPECT_THAT(RelativeError(keplerian_elements.inclination,
                            keplerian_elements_1960_.inclination),
              IsNear(9.3e-10));
  EXPECT_THAT(
      RelativeError(keplerian_elements.longitude_of_ascending_node,
                    keplerian_elements_1960_.longitude_of_ascending_node),
      IsNear(7.5e-9));
  EXPECT_THAT(*keplerian_elements_1960_.argument_of_periapsis -
                  *keplerian_elements.argument_of_periapsis,
              IsNear(4.206 * ArcSecond, 1.0002));
  EXPECT_THAT(*keplerian_elements.mean_anomaly -
                  *keplerian_elements_1960_.mean_anomaly,
              IsNear(17.03 * ArcSecond, 1.0003));
}

TEST_F(MercuryPerihelionTest, DISABLED_Year2050) {
  ephemeris_->Prolong(t_2050_);

  auto const& sun_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Sun");
  auto const& mercury_trajectory =
      solar_system_1950_.trajectory(*ephemeris_, "Mercury");

  RelativeDegreesOfFreedom<ICRS> const relative_degrees_of_freedom =
      mercury_trajectory.EvaluateDegreesOfFreedom(t_2050_) -
      sun_trajectory.EvaluateDegreesOfFreedom(t_2050_);
  KeplerOrbit<ICRS> orbit(
      *sun_, *mercury_, relative_degrees_of_freedom, t_2050_);
  KeplerianElements<ICRS> const keplerian_elements = orbit.elements_at_epoch();

  EXPECT_THAT(RelativeError(*keplerian_elements.eccentricity,
                            *keplerian_elements_2050_.eccentricity),
              IsNear(9.8e-8));
  EXPECT_THAT(RelativeError(*keplerian_elements.semimajor_axis,
                            *keplerian_elements_2050_.semimajor_axis),
              IsNear(1.8e-8));
  EXPECT_THAT(RelativeError(*keplerian_elements.mean_motion,
                            *keplerian_elements_2050_.mean_motion),
              IsNear(2.7e-8));
  EXPECT_THAT(RelativeError(keplerian_elements.inclination,
                            keplerian_elements_2050_.inclination),
              IsNear(4.0e-9));
  EXPECT_THAT(
      RelativeError(keplerian_elements.longitude_of_ascending_node,
                    keplerian_elements_2050_.longitude_of_ascending_node),
      IsNear(1.5e-7));
  // The actual number is 42.98" on average, but there are other periodic
  // effects so your mileage will vary.  See Nobili and Will, The real value of
  // Mercury's perihelion advance.
  EXPECT_THAT(*keplerian_elements_2050_.argument_of_periapsis -
                  *keplerian_elements.argument_of_periapsis,
              IsNear(42.87 * ArcSecond, 1.0001));
  EXPECT_THAT(*keplerian_elements.mean_anomaly -
                  *keplerian_elements_2050_.mean_anomaly,
              IsNear(171.2 * ArcSecond, 1.0006));
}

#endif

}  // namespace astronomy
}  // namespace principia
