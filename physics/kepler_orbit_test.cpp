
#include "physics/kepler_orbit.hpp"

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/time_scales.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace physics {
namespace internal_kepler_orbit {

using astronomy::ICRFJ2000Equator;
using astronomy::J2000;
using astronomy::operator""_TT;
using quantities::astronomy::JulianYear;
using quantities::si::AstronomicalUnit;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;

// Target body name : Moon(301) { source: DE431mx }
// Centre body name : Earth(399) { source: DE431mx }
// Output units    : KM-S, deg, Julian day number (Tp)
// Reference frame : ICRF / J2000.0
// Coordinate systm : Earth Mean Equator and Equinox of Reference Epoch
// System GM : 4.0350323550225975e+05 km^3/s^2
// 2457397.500000000 = A.D. 2016-Jan-10 00:00:00.0000 (TDB)
// EC= 4.772161502830355e-02 QR= 3.685366825859605e+05 IN= 1.842335956339145e+01
// OM= 1.752118723367974e+00 W = 3.551364385683149e+02 Tp=  2457402.376861531753
// N = 1.511718576836574e-04 MA= 2.963020996150547e+02 TA= 2.912732951134421e+02
// A = 3.870051955415476e+05 AD= 4.054737084971346e+05 PR= 2.381395621619845e+06
// X = 1.177367562036580e+05 Y =-3.419908628150604e+05 Z =-1.150659799281941e+05
// VX= 9.745048087261129e-01 VY= 3.500672337210811e-01 VZ= 1.066306010215636e-01
// Symbol meaning
//   EC     Eccentricity, e
//   QR     Periapsis distance, q (km)
//   IN     Inclination w.r.t xy-plane, i (degrees)
//   OM     Longitude of Ascending Node, OMEGA, (degrees)
//   W      Argument of Perifocus, w (degrees)
//   Tp     Time of periapsis (Julian day number)
//   N      Mean motion, n (degrees/sec)
//   MA     Mean anomaly, M (degrees)
//   TA     True anomaly, nu (degrees)
//   A      Semi-major axis, a (km)
//   AD     Apoapsis distance (km)
//   PR     Sidereal orbit period (sec)

// Target body name: Voyager 1 (spacecraft) (-31)    {source: Voyager_1}
// Center body name: Sun (10)                        {source: Voyager_1}
// 2457899.322222222 = A.D. 2017-May-25 19:44:00.0000 TDB
// EC= 3.754904752975423e+00 QR= 1.324687933572433e+09 IN= 1.231535474702641e+01
// OM= 1.762814430848047e+02 W = 3.412521452280985e+02 Tp=  2444230.840535562020
// N = 1.979556771581467e-06 MA= 2.337771065477673e+03 TA= 1.006849541017370e+02
// A =-4.808470899553643e+08 AD= 9.999999000000000e+99 PR= 9.999999000000000e+99
// X =-4.202547896125371e+09 Y =-1.982453908150731e+10 Z = 4.378406169173994e+09
// VX=-2.067668772297011e+00 VY=-1.647515322877371e+01 VZ= 3.618493348898172e+00


class KeplerOrbitTest : public ::testing::Test {
 protected:
  static KeplerianElements<ICRFJ2000Equator> MoonElements() {
    KeplerianElements<ICRFJ2000Equator> elements;
    elements.eccentricity                = 4.772161502830355e-02;
    elements.semimajor_axis              = 3.870051955415476e+05 * Kilo(Metre);
    elements.mean_motion                 = 1.511718576836574e-04 * (Degree /
                                                                    Second);
    elements.period                      = 2.381395621619845e+06 * Second;
    elements.periapsis_distance          = 3.685366825859605e+05 * Kilo(Metre);
    elements.apoapsis_distance           = 4.054737084971346e+05 * Kilo(Metre);
    elements.inclination                 = 1.842335956339145e+01 * Degree;
    elements.longitude_of_ascending_node = 1.752118723367974e+00 * Degree;
    elements.argument_of_periapsis       = 3.551364385683149e+02 * Degree;
    elements.mean_anomaly                = 2.963020996150547e+02 * Degree;
    elements.true_anomaly                = 2.912732951134421e+02 * Degree;
    return elements;
  }

  static KeplerianElements<ICRFJ2000Equator> VoyagerElements() {
    KeplerianElements<ICRFJ2000Equator> elements;
    elements.eccentricity                =  3.754904752975423e+00;
    elements.semimajor_axis              = -4.808470899553643e+08 * Kilo(Metre);
    elements.hyperbolic_mean_motion      =  1.979556771581467e-06 * (Degree /
                                                                     Second);
    elements.periapsis_distance          =  1.324687933572433e+09 * Kilo(Metre);
    elements.inclination                 =  1.231535474702641e+01 * Degree;
    elements.longitude_of_ascending_node =  1.762814430848047e+02 * Degree;
    elements.argument_of_periapsis       =  3.412521452280985e+02 * Degree;
    elements.hyperbolic_mean_anomaly     =  2.337771065477673e+03 * Degree;
    elements.true_anomaly                =  1.006849541017370e+02 * Degree;
    return elements;
  }

  static KeplerianElements<ICRFJ2000Equator> SimpleEllipse() {
    KeplerianElements<ICRFJ2000Equator> elements;
    elements.eccentricity = 0.5;
    elements.asymptotic_true_anomaly = -NaN<Angle>();
    elements.turning_angle = -NaN<Angle>();
    elements.semimajor_axis = 1 * AstronomicalUnit;
    elements.specific_energy =
        -0.5 * Pow<2>(AstronomicalUnit) / Pow<2>(JulianYear);
    elements.characteristic_energy =
        -1 * Pow<2>(AstronomicalUnit) / Pow<2>(JulianYear);
    elements.period = 2 * π * JulianYear;
    elements.mean_motion = 1 * Radian / JulianYear;
    elements.hyperbolic_mean_motion = -NaN<AngularFrequency>();
    elements.hyperbolic_excess_velocity = -NaN<Speed>();
    elements.semiminor_axis = Sqrt(3) / 2 * AstronomicalUnit;
    elements.impact_parameter = -NaN<Length>();
    elements.semilatus_rectum = 0.75 * AstronomicalUnit;
    elements.specific_angular_momentum =
        (Sqrt(3) / 2 * Pow<2>(AstronomicalUnit) / JulianYear) * Radian;
    elements.periapsis_distance = 0.5 * AstronomicalUnit;
    elements.apoapsis_distance = 1.5 * AstronomicalUnit;

    elements.inclination = 1 * Degree;
    elements.longitude_of_ascending_node = 0.5 * Radian;
    elements.argument_of_periapsis = 0.5 * Radian;
    elements.longitude_of_periapsis = 1 * Radian;

    elements.true_anomaly = π / 2 * Radian;
    elements.mean_anomaly = (π / 3 - Sqrt(3) / 4) * Radian;
    elements.hyperbolic_mean_anomaly = +NaN<Angle>();
    return elements;
  }

  static KeplerianElements<ICRFJ2000Equator> SimpleHyperbola() {
    KeplerianElements<ICRFJ2000Equator> elements;
    elements.eccentricity = 1.5;
    elements.asymptotic_true_anomaly = ArcCos(-1 / 1.5);
    elements.turning_angle = 2 * ArcSin(1 / 1.5);
    elements.semimajor_axis = -1 * AstronomicalUnit;
    elements.specific_energy =
        0.5 * Pow<2>(AstronomicalUnit) / Pow<2>(JulianYear);
    elements.characteristic_energy =
        1 * Pow<2>(AstronomicalUnit) / Pow<2>(JulianYear);
    elements.period = -NaN<Time>();
    elements.mean_motion = -NaN<AngularFrequency>();
    elements.hyperbolic_mean_motion = 1 * Radian / JulianYear;
    elements.hyperbolic_excess_velocity = 1 * AstronomicalUnit / JulianYear;
    elements.semiminor_axis = -NaN<Length>();
    elements.impact_parameter = Sqrt(5) / 2 * AstronomicalUnit;
    elements.semilatus_rectum = 1.25 * AstronomicalUnit;
    elements.specific_angular_momentum =
        (Sqrt(5) / 2 * Pow<2>(AstronomicalUnit) / JulianYear) * Radian;
    elements.periapsis_distance = 0.5 * AstronomicalUnit;
    elements.apoapsis_distance = -2.5 * AstronomicalUnit;

    elements.inclination = 1 * Degree;
    elements.longitude_of_ascending_node = 0.5 * Radian;
    elements.argument_of_periapsis = 0.5 * Radian;
    elements.longitude_of_periapsis = 1 * Radian;

    elements.true_anomaly = π / 2 * Radian;
    elements.mean_anomaly = -NaN<Angle>();
    elements.hyperbolic_mean_anomaly =
        (3 * Sqrt(5) / 4) * Radian - ArcCosh(1.5);
    return elements;
  }

  // An ellipse with a = 1 au, e very close to 1.
  static KeplerianElements<ICRFJ2000Equator> NearlyParabolicEllipse() {
    KeplerianElements<ICRFJ2000Equator> elements;
    // This is about half the bits at 1, which maximizes the error in the naïve
    // evaluation of 1 - (1 - ε)².  Note the 1 - (1 - x), to ensure that
    // (1 - ε) is exact.
    // Note that [1 - (1 - ε)²] = [[2ε] - [ε²]], where [] denotes rounding.
    constexpr double ε = 1 - (1 - 2e-9);
    constexpr double ε² = ε * ε;
    constexpr double ε³ = ε² * ε;
    constexpr double ε⁴ = ε² * ε²;
    elements.eccentricity = 1 - ε;
    elements.asymptotic_true_anomaly = -NaN<Angle>();
    elements.turning_angle = -NaN<Angle>();
    elements.semimajor_axis = 1 * AstronomicalUnit;
    elements.specific_energy =
        -0.5 * Pow<2>(AstronomicalUnit) / Pow<2>(JulianYear);
    elements.characteristic_energy =
        -1 * Pow<2>(AstronomicalUnit) / Pow<2>(JulianYear);
    elements.period = 2 * π * JulianYear;
    elements.mean_motion = 1 * Radian / JulianYear;
    elements.hyperbolic_mean_motion = -NaN<AngularFrequency>();
    elements.hyperbolic_excess_velocity = -NaN<Speed>();
    elements.semiminor_axis = Sqrt(2 * ε - ε²) * AstronomicalUnit;
    elements.impact_parameter = -NaN<Length>();
    elements.semilatus_rectum = (2 * ε - ε²) * AstronomicalUnit;
    elements.specific_angular_momentum =
        (Sqrt(2 * ε - ε²) * Pow<2>(AstronomicalUnit) / JulianYear) * Radian;
    elements.periapsis_distance = ε * AstronomicalUnit;
    elements.apoapsis_distance = (2 - ε) * AstronomicalUnit;

    elements.inclination = 1 * Degree;
    elements.longitude_of_ascending_node = 0.5 * Radian;
    elements.argument_of_periapsis = 0.5 * Radian;
    elements.longitude_of_periapsis = 1 * Radian;

    elements.true_anomaly = π / 2 * Radian;
    // This expression gives the correctly-rounded result for a true anomaly of
    // π / 2, and is off by 0.96 ULPs for a true anomaly of [π / 2].
    elements.mean_anomaly = Sqrt(32 * ε³ / 9 - 16 * ε⁴ / 15) * Radian;
    elements.hyperbolic_mean_anomaly = +NaN<Angle>();
    return elements;
  }

  static void ExpectConicParametersAlmostEqual(
      KeplerianElements<ICRFJ2000Equator> const& actual,
      KeplerianElements<ICRFJ2000Equator> const& expected,
      std::int64_t const eccentrity_ulps,
      std::int64_t const asymptotic_true_anomaly_ulps,
      std::int64_t const turning_angle_ulps,
      std::int64_t const semimajor_axis_ulps,
      std::int64_t const specific_energy_ulps,
      std::int64_t const characteristic_energy_ulps,
      std::int64_t const mean_motion_ulps,
      std::int64_t const period_ulps,
      std::int64_t const hyperbolic_mean_motion_ulps,
      std::int64_t const hyperbolic_excess_velocity_ulps,
      std::int64_t const semiminor_axis_ulps,
      std::int64_t const impact_parameter_ulps,
      std::int64_t const semilatus_rectum_ulps,
      std::int64_t const specific_angular_momentum_ulps,
      std::int64_t const periapsis_distance_ulps,
      std::int64_t const apoapsis_distance_ulps) {
    EXPECT_THAT(*actual.eccentricity,
                AlmostEquals(*expected.eccentricity, eccentrity_ulps));
    EXPECT_THAT(*actual.asymptotic_true_anomaly,
                AlmostEquals(*expected.asymptotic_true_anomaly,
                             asymptotic_true_anomaly_ulps));
    EXPECT_THAT(*actual.turning_angle,
                AlmostEquals(*expected.turning_angle, turning_angle_ulps));
    EXPECT_THAT(*actual.semimajor_axis,
                AlmostEquals(*expected.semimajor_axis, semimajor_axis_ulps));
    EXPECT_THAT(*actual.specific_energy,
                AlmostEquals(*expected.specific_energy, specific_energy_ulps));
    EXPECT_THAT(*actual.characteristic_energy,
                AlmostEquals(*expected.characteristic_energy,
                             characteristic_energy_ulps));
    EXPECT_THAT(*actual.mean_motion,
                AlmostEquals(*expected.mean_motion, mean_motion_ulps));
    EXPECT_THAT(*actual.period, AlmostEquals(*expected.period, period_ulps));
    EXPECT_THAT(*actual.hyperbolic_mean_motion,
                AlmostEquals(*expected.hyperbolic_mean_motion,
                             hyperbolic_mean_motion_ulps));
    EXPECT_THAT(*actual.hyperbolic_excess_velocity,
                AlmostEquals(*expected.hyperbolic_excess_velocity,
                             hyperbolic_excess_velocity_ulps));
    EXPECT_THAT(*actual.semiminor_axis,
                AlmostEquals(*expected.semiminor_axis, semiminor_axis_ulps));
    EXPECT_THAT(
        *actual.impact_parameter,
        AlmostEquals(*expected.impact_parameter, impact_parameter_ulps));
    EXPECT_THAT(
        *actual.semilatus_rectum,
        AlmostEquals(*expected.semilatus_rectum, semilatus_rectum_ulps));
    EXPECT_THAT(*actual.specific_angular_momentum,
                AlmostEquals(*expected.specific_angular_momentum,
                             specific_angular_momentum_ulps));
    EXPECT_THAT(
        *actual.periapsis_distance,
        AlmostEquals(*expected.periapsis_distance, periapsis_distance_ulps));
    EXPECT_THAT(
        *actual.apoapsis_distance,
        AlmostEquals(*expected.apoapsis_distance, apoapsis_distance_ulps));
  }

  MassiveBody const body_{1 * Pow<3>(AstronomicalUnit) / Pow<2>(JulianYear)};
};

TEST_F(KeplerOrbitTest, EarthMoon) {
  SolarSystem<ICRFJ2000Equator> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2433282_500000000.proto.txt");
  auto const earth = SolarSystem<ICRFJ2000Equator>::MakeMassiveBody(
                         solar_system.gravity_model_message("Earth"));
  auto const moon = SolarSystem<ICRFJ2000Equator>::MakeMassiveBody(
                        solar_system.gravity_model_message("Moon"));
  // The numbers in the gravity models and those from the query above both come
  // from DE431, so the sums are the same up to round-off.
  EXPECT_THAT(
      earth->gravitational_parameter() + moon->gravitational_parameter(),
      AlmostEquals(
          4.0350323550225975e+05 * (Pow<3>(Kilo(Metre)) / Pow<2>(Second)), 1));
  constexpr Instant date = "JD2457397.500000000"_TT;

  Displacement<ICRFJ2000Equator> const expected_displacement(
      { 1.177367562036580e+05 * Kilo(Metre),
       -3.419908628150604e+05 * Kilo(Metre),
       -1.150659799281941e+05 * Kilo(Metre)});
  Velocity<ICRFJ2000Equator> const expected_velocity(
      {9.745048087261129e-01 * (Kilo(Metre) / Second),
       3.500672337210811e-01 * (Kilo(Metre) / Second),
       1.066306010215636e-01 * (Kilo(Metre) / Second)});

  auto partial_elements = MoonElements();
  partial_elements.mean_motion.reset();
  partial_elements.period.reset();
  partial_elements.periapsis_distance.reset();
  partial_elements.apoapsis_distance.reset();
  partial_elements.true_anomaly.reset();
  {
    KeplerOrbit<ICRFJ2000Equator> moon_orbit(
        *earth, *moon, partial_elements, date);
    EXPECT_THAT(moon_orbit.StateVectors(date).displacement(),
                AlmostEquals(expected_displacement, 15));
    EXPECT_THAT(moon_orbit.StateVectors(date).velocity(),
                AlmostEquals(expected_velocity, 21));
    EXPECT_THAT(*moon_orbit.elements_at_epoch().mean_motion,
                AlmostEquals(*MoonElements().mean_motion, 2));
  }

  partial_elements.semimajor_axis.reset();
  partial_elements.periapsis_distance = MoonElements().periapsis_distance;
  {
    KeplerOrbit<ICRFJ2000Equator> moon_orbit(
        *earth, *moon, partial_elements, date);
    EXPECT_THAT(moon_orbit.StateVectors(date).displacement(),
                AlmostEquals(expected_displacement, 13, 15));
    EXPECT_THAT(moon_orbit.StateVectors(date).velocity(),
                AlmostEquals(expected_velocity, 23));
  }

  KeplerOrbit<ICRFJ2000Equator> moon_orbit(
      *earth, *moon, {expected_displacement, expected_velocity}, date);
  EXPECT_THAT(*moon_orbit.elements_at_epoch().eccentricity,
              AlmostEquals(*MoonElements().eccentricity, 8));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().semimajor_axis,
              AlmostEquals(*MoonElements().semimajor_axis, 1));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().mean_motion,
              AlmostEquals(*MoonElements().mean_motion, 0));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().period,
              AlmostEquals(*MoonElements().period, 0));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().periapsis_distance,
              AlmostEquals(*MoonElements().periapsis_distance, 2));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().apoapsis_distance,
              AlmostEquals(*MoonElements().apoapsis_distance, 1));
  EXPECT_THAT(moon_orbit.elements_at_epoch().inclination,
              AlmostEquals(MoonElements().inclination, 1));
  EXPECT_THAT(moon_orbit.elements_at_epoch().longitude_of_ascending_node,
              AlmostEquals(MoonElements().longitude_of_ascending_node, 28));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().argument_of_periapsis,
              AlmostEquals(*MoonElements().argument_of_periapsis, 6));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().mean_anomaly,
              AlmostEquals(*MoonElements().mean_anomaly, 6));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().true_anomaly,
              AlmostEquals(*MoonElements().true_anomaly, 5));
}

TEST_F(KeplerOrbitTest, Voyager1) {
  SolarSystem<ICRFJ2000Equator> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2433282_500000000.proto.txt");
  auto const sun = SolarSystem<ICRFJ2000Equator>::MakeMassiveBody(
                         solar_system.gravity_model_message("Sun"));
  MasslessBody const voyager1{};
  constexpr Instant date = "2017-05-25T19:44:00,000"_TT;

  Displacement<ICRFJ2000Equator> const expected_displacement(
      {-4.202547896125371e+09 * Kilo(Metre),
       -1.982453908150731e+10 * Kilo(Metre),
        4.378406169173994e+09 * Kilo(Metre)});
  Velocity<ICRFJ2000Equator> const expected_velocity(
      {-2.067668772297011e+00 * (Kilo(Metre) / Second),
       -1.647515322877371e+01 * (Kilo(Metre) / Second),
        3.618493348898172e+00 * (Kilo(Metre) / Second)});

  auto partial_elements = VoyagerElements();
  partial_elements.hyperbolic_mean_motion.reset();
  partial_elements.periapsis_distance.reset();
  partial_elements.true_anomaly.reset();
  {
    KeplerOrbit<ICRFJ2000Equator> voyager_orbit(
        *sun, voyager1, partial_elements, date);
    EXPECT_THAT(voyager_orbit.StateVectors(date).displacement(),
                AlmostEquals(expected_displacement, 37, 49));
    EXPECT_THAT(voyager_orbit.StateVectors(date).velocity(),
                AlmostEquals(expected_velocity, 26));
    EXPECT_THAT(*voyager_orbit.elements_at_epoch().hyperbolic_mean_motion,
                AlmostEquals(*VoyagerElements().hyperbolic_mean_motion, 4));
  }

  partial_elements.semimajor_axis.reset();
  partial_elements.periapsis_distance = VoyagerElements().periapsis_distance;
  {
    KeplerOrbit<ICRFJ2000Equator> voyager_orbit(
        *sun, voyager1, partial_elements, date);
    EXPECT_THAT(voyager_orbit.StateVectors(date).displacement(),
                AlmostEquals(expected_displacement, 31, 44));
    EXPECT_THAT(voyager_orbit.StateVectors(date).velocity(),
                AlmostEquals(expected_velocity, 27, 28));
  }

  KeplerOrbit<ICRFJ2000Equator> voyager_orbit(
      *sun, voyager1, {expected_displacement, expected_velocity}, date);
  EXPECT_THAT(*voyager_orbit.elements_at_epoch().eccentricity,
              AlmostEquals(*VoyagerElements().eccentricity, 2));
  EXPECT_THAT(*voyager_orbit.elements_at_epoch().semimajor_axis,
              AlmostEquals(*VoyagerElements().semimajor_axis, 3));
  EXPECT_THAT(*voyager_orbit.elements_at_epoch().hyperbolic_mean_motion,
              AlmostEquals(*VoyagerElements().hyperbolic_mean_motion, 0));
  EXPECT_THAT(*voyager_orbit.elements_at_epoch().periapsis_distance,
              AlmostEquals(*VoyagerElements().periapsis_distance, 4));
  EXPECT_THAT(voyager_orbit.elements_at_epoch().inclination,
              AlmostEquals(VoyagerElements().inclination, 10));
  EXPECT_THAT(voyager_orbit.elements_at_epoch().longitude_of_ascending_node,
              AlmostEquals(VoyagerElements().longitude_of_ascending_node, 9));
  EXPECT_THAT(*voyager_orbit.elements_at_epoch().argument_of_periapsis,
              AlmostEquals(*VoyagerElements().argument_of_periapsis, 5));
  EXPECT_THAT(*voyager_orbit.elements_at_epoch().hyperbolic_mean_anomaly,
              AlmostEquals(*VoyagerElements().hyperbolic_mean_anomaly, 0));
  EXPECT_THAT(*voyager_orbit.elements_at_epoch().true_anomaly,
              AlmostEquals(*VoyagerElements().true_anomaly, 3));
}

TEST_F(KeplerOrbitTest, TrueAnomalyToEllipticMeanAnomaly) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semilatus_rectum = SimpleEllipse().semilatus_rectum;
  elements.periapsis_distance = SimpleEllipse().periapsis_distance;
  elements.argument_of_periapsis.emplace();
  elements.true_anomaly = SimpleEllipse().true_anomaly;
  elements =
      KeplerOrbit<ICRFJ2000Equator>(body_, MasslessBody{}, elements, J2000)
          .elements_at_epoch();
  EXPECT_THAT(*elements.true_anomaly, Eq(*SimpleEllipse().true_anomaly));
  EXPECT_THAT(*elements.mean_anomaly,
              AlmostEquals(*SimpleEllipse().mean_anomaly, 0));
  EXPECT_THAT(*elements.hyperbolic_mean_anomaly,
              AlmostEquals(*SimpleEllipse().hyperbolic_mean_anomaly, 0));
}

TEST_F(KeplerOrbitTest, TrueAnomalyToHyperbolicMeanAnomaly) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semilatus_rectum = SimpleHyperbola().semilatus_rectum;
  elements.periapsis_distance = SimpleHyperbola().periapsis_distance;
  elements.argument_of_periapsis.emplace();
  elements.true_anomaly = SimpleHyperbola().true_anomaly;
  elements =
      KeplerOrbit<ICRFJ2000Equator>(body_, MasslessBody{}, elements, J2000)
          .elements_at_epoch();
  EXPECT_THAT(*elements.true_anomaly, Eq(*SimpleHyperbola().true_anomaly));
  EXPECT_THAT(*elements.mean_anomaly,
              AlmostEquals(*SimpleHyperbola().mean_anomaly, 0));
  EXPECT_THAT(*elements.hyperbolic_mean_anomaly,
              AlmostEquals(*SimpleHyperbola().hyperbolic_mean_anomaly, 0));
}

TEST_F(KeplerOrbitTest, EllipticMeanAnomalyToTrueAnomaly) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semilatus_rectum = SimpleEllipse().semilatus_rectum;
  elements.periapsis_distance = SimpleEllipse().periapsis_distance;
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly = SimpleEllipse().mean_anomaly;
  elements =
      KeplerOrbit<ICRFJ2000Equator>(body_, MasslessBody{}, elements, J2000)
          .elements_at_epoch();
  EXPECT_THAT(*elements.mean_anomaly, Eq(*SimpleEllipse().mean_anomaly));
  EXPECT_THAT(*elements.true_anomaly,
              AlmostEquals(*SimpleEllipse().true_anomaly, 1));
  EXPECT_THAT(*elements.hyperbolic_mean_anomaly,
              AlmostEquals(*SimpleEllipse().hyperbolic_mean_anomaly, 0));
}

TEST_F(KeplerOrbitTest, HyperbolicMeanAnomalyToTrueAnomaly) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semilatus_rectum = SimpleHyperbola().semilatus_rectum;
  elements.periapsis_distance = SimpleHyperbola().periapsis_distance;
  elements.argument_of_periapsis.emplace();
  elements.hyperbolic_mean_anomaly = SimpleHyperbola().hyperbolic_mean_anomaly;
  elements =
      KeplerOrbit<ICRFJ2000Equator>(body_, MasslessBody{}, elements, J2000)
          .elements_at_epoch();
  EXPECT_THAT(*elements.hyperbolic_mean_anomaly,
              Eq(*SimpleHyperbola().hyperbolic_mean_anomaly));
  EXPECT_THAT(*elements.true_anomaly,
              AlmostEquals(*SimpleHyperbola().true_anomaly, 0));
  EXPECT_THAT(*elements.mean_anomaly,
              AlmostEquals(*SimpleHyperbola().mean_anomaly, 0));
}

TEST_F(KeplerOrbitTest, OrientationFromLongitudeOfPeriapsis) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity = 0;
  elements.semimajor_axis = 1 * Metre;
  elements.mean_anomaly.emplace();
  elements.inclination = SimpleEllipse().inclination;
  elements.longitude_of_ascending_node =
      SimpleEllipse().longitude_of_ascending_node;
  elements.longitude_of_periapsis = SimpleEllipse().longitude_of_periapsis;
  elements =
      KeplerOrbit<ICRFJ2000Equator>(body_, MasslessBody{}, elements, J2000)
          .elements_at_epoch();
  EXPECT_THAT(*elements.longitude_of_periapsis,
              Eq(*SimpleEllipse().longitude_of_periapsis));
  EXPECT_THAT(*elements.argument_of_periapsis,
              AlmostEquals(*SimpleEllipse().argument_of_periapsis, 0));
}

TEST_F(KeplerOrbitTest, OrientationFromArgumentOfPeriapsis) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity = 0;
  elements.semimajor_axis = 1 * Metre;
  elements.mean_anomaly.emplace();
  elements.inclination = SimpleEllipse().inclination;
  elements.longitude_of_ascending_node =
      SimpleEllipse().longitude_of_ascending_node;
  elements.argument_of_periapsis = SimpleEllipse().argument_of_periapsis;
  elements =
      KeplerOrbit<ICRFJ2000Equator>(body_, MasslessBody{}, elements, J2000)
          .elements_at_epoch();
  EXPECT_THAT(*elements.argument_of_periapsis,
              Eq(*SimpleEllipse().argument_of_periapsis));
  EXPECT_THAT(*elements.longitude_of_periapsis,
              AlmostEquals(*SimpleEllipse().longitude_of_periapsis, 0));
}

#define CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(element1, element2, reference)       \
  \
[&]() {                                                                        \
    KeplerianElements<ICRFJ2000Equator> elements;                              \
    elements.element1 = (reference).element1;                                  \
    elements.element2 = (reference).element2;                                  \
    /* Leaving the orientation parameters and anomalies to their default    */ \
    /* values.  This test does not exercise them.                           */ \
    elements.argument_of_periapsis.emplace();                                  \
    elements.mean_anomaly.emplace();                                           \
    elements =                                                                 \
        KeplerOrbit<ICRFJ2000Equator>(body_, MasslessBody{}, elements, J2000)  \
            .elements_at_epoch();                                              \
    /* The inputs must not change.                                          */ \
    EXPECT_THAT(*elements.element1, Eq(*(reference).element1));                \
    EXPECT_THAT(*elements.element2, Eq(*(reference).element2));                \
    return elements;                                                           \
  }()

// Test all code paths for completing the conic parameters, first with an
// ellipse.

// Test all choices of two categories of conic parameters.

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndSemimajorAxis) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, semimajor_axis, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndSemiminorAxis) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, semiminor_axis, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/2,
                                   /*characteristic_energy_ulps=*/2,
                                   /*mean_motion_ulps=*/3,
                                   /*period_ulps=*/3,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/2);
}

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, semilatus_rectum, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, periapsis_distance, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, apoapsis_distance, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromSemimajorAxisAndSemiminorAxis) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semimajor_axis, semiminor_axis, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/2,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/2,
                                   /*apoapsis_distance_ulps=*/1);
}

TEST_F(KeplerOrbitTest, EllipseFromSemimajorAxisAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semimajor_axis, semilatus_rectum, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromSemimajorAxisAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semimajor_axis, periapsis_distance, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromSemimajorAxisAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semimajor_axis, apoapsis_distance, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromSemiminorAxisAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semiminor_axis, semilatus_rectum, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/2,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/2,
                                   /*characteristic_energy_ulps=*/2,
                                   /*mean_motion_ulps=*/3,
                                   /*period_ulps=*/3,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/2);
}

TEST_F(KeplerOrbitTest, EllipseFromSemiminorAxisAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semiminor_axis, periapsis_distance, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/2,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/1);
}

TEST_F(KeplerOrbitTest, EllipseFromSemiminorAxisAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semiminor_axis, apoapsis_distance, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/1,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromSemilatusRectumAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semilatus_rectum, periapsis_distance, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromSemilatusRectumAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semilatus_rectum, apoapsis_distance, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, EllipseFromPeriapsisDistanceAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          periapsis_distance, apoapsis_distance, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

// Test all alternative semimajor axis specifications (the semimajor axis itself
// is already tested above).

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndSpecificEnergy) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, specific_energy, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/0,
                                   /*characteristic_energy_ulps=*/0,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/1,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/2,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/2,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/2);
}

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndCharacteristicEnergy) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, characteristic_energy, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/0,
                                   /*characteristic_energy_ulps=*/0,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/1,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/2,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/2,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/2);
}

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndMeanMotion) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, mean_motion, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/2,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/2,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/2);
}

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndPeriod) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(eccentricity, period, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/0,
                                   /*characteristic_energy_ulps=*/0,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

// Test all alternative semilatus rectum specifications (the semilatus rectum
// itself is already tested above).

TEST_F(KeplerOrbitTest, EllipseFromEccentricityAndSpecificAngularMomentum) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, specific_angular_momentum, SimpleEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

// Similar tests with a hyperbola.

// Test all choices of two categories of conic parameters.

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndSemimajorAxis) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, semimajor_axis, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndImpactParameter) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, impact_parameter, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/0,
                                   /*characteristic_energy_ulps=*/0,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/1);
}

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, semilatus_rectum, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, periapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, apoapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromSemimajorAxisAndImpactParameter) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semimajor_axis, impact_parameter, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromSemimajorAxisAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semimajor_axis, semilatus_rectum, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromSemimajorAxisAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semimajor_axis, periapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromSemimajorAxisAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semimajor_axis, apoapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromImpactParameterAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          impact_parameter, semilatus_rectum, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromImpactParameterAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          impact_parameter, periapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/1);
}

TEST_F(KeplerOrbitTest, HyperbolaFromImpactParameterAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          impact_parameter, apoapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromSemilatusRectumAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semilatus_rectum, periapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromSemilatusRectumAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          semilatus_rectum, apoapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest, HyperbolaFromPeriapsisDistanceAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          periapsis_distance, apoapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

// Test all alternative eccentricity specifications (the eccentricity itself is
// already tested above).

TEST_F(KeplerOrbitTest,
       HyperbolaFromAsymptoticTrueAnomalyAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          asymptotic_true_anomaly, periapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/1,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/2,
                                   /*semimajor_axis_ulps=*/2,
                                   /*specific_energy_ulps=*/2,
                                   /*characteristic_energy_ulps=*/2,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/2,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/2);
}

TEST_F(KeplerOrbitTest, HyperbolaFromTurningAngleAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          turning_angle, periapsis_distance, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/1,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

// Test all alternative semimajor axis specifications (the semimajor axis itself
// is already tested above).

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndSpecificEnergy) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, specific_energy, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/0,
                                   /*characteristic_energy_ulps=*/0,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/1);
}

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndCharacteristicEnergy) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, characteristic_energy, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/0,
                                   /*characteristic_energy_ulps=*/0,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/1);
}

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndHyperbolicMeanMotion) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, hyperbolic_mean_motion, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
#if PRINCIPIA_COMPILER_MSVC
                                   /*hyperbolic_excess_velocity_ulps=*/1,
#else
                                   /*hyperbolic_excess_velocity_ulps=*/0,
#endif
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/1);
}

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndHyperbolicExcessVelocity) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, hyperbolic_excess_velocity, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/2,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/2,
                                   /*semilatus_rectum_ulps=*/2,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/2,
                                   /*apoapsis_distance_ulps=*/2);
}

// Test all alternative semilatus rectum specifications (the semilatus rectum
// itself is already tested above).

TEST_F(KeplerOrbitTest, HyperbolaFromEccentricityAndSpecificAngularMomentum) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, specific_angular_momentum, SimpleHyperbola());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/SimpleHyperbola(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/1,
                                   /*specific_energy_ulps=*/0,
                                   /*characteristic_energy_ulps=*/0,
                                   /*mean_motion_ulps=*/0,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/1,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/1,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/1,
                                   /*apoapsis_distance_ulps=*/1);
}

// Tests that emphasize ill-conditioning.

TEST_F(KeplerOrbitTest,
       NearlyParabolicEllipseFromEccentricityAndSemimajorAxis) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, semimajor_axis, NearlyParabolicEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/NearlyParabolicEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/2539776,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/5263508,
                                   /*specific_angular_momentum_ulps=*/2939388,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest,
       NearlyParabolicEllipseFromEccentricityAndSemiminorAxis) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, semiminor_axis, NearlyParabolicEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/NearlyParabolicEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/2451012,
                                   /*specific_energy_ulps=*/3016150,
                                   /*characteristic_energy_ulps=*/3016150,
                                   /*mean_motion_ulps=*/3591428,
                                   /*period_ulps=*/4989938,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/0,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/2631754,
                                   /*specific_angular_momentum_ulps=*/1469694,
                                   /*periapsis_distance_ulps=*/2631754,
                                   /*apoapsis_distance_ulps=*/2451012);
}

TEST_F(KeplerOrbitTest,
       NearlyParabolicEllipseFromEccentricityAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, semilatus_rectum, NearlyParabolicEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/NearlyParabolicEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/4902024,
                                   /*specific_energy_ulps=*/6032300,
                                   /*characteristic_energy_ulps=*/6032300,
                                   /*mean_motion_ulps=*/7182855,
                                   /*period_ulps=*/9979876,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/2539775,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/1);
}

TEST_F(KeplerOrbitTest,
       NearlyParabolicEllipseFromEccentricityAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, periapsis_distance, NearlyParabolicEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/NearlyParabolicEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/2539776,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/0,
                                   /*specific_angular_momentum_ulps=*/0,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

TEST_F(KeplerOrbitTest,
       NearlyParabolicEllipseFromEccentricityAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> const elements =
      CONSTRUCT_CONIC_FROM_TWO_ELEMENTS(
          eccentricity, apoapsis_distance, NearlyParabolicEllipse());
  ExpectConicParametersAlmostEqual(/*actual=*/elements,
                                   /*expected=*/NearlyParabolicEllipse(),
                                   /*eccentrity_ulps=*/0,
                                   /*asymptotic_true_anomaly_ulps=*/0,
                                   /*turning_angle_ulps=*/0,
                                   /*semimajor_axis_ulps=*/0,
                                   /*specific_energy_ulps=*/1,
                                   /*characteristic_energy_ulps=*/1,
                                   /*mean_motion_ulps=*/1,
                                   /*period_ulps=*/0,
                                   /*hyperbolic_mean_motion_ulps=*/0,
                                   /*hyperbolic_excess_velocity_ulps=*/0,
                                   /*semiminor_axis_ulps=*/2539776,
                                   /*impact_parameter_ulps=*/0,
                                   /*semilatus_rectum_ulps=*/1,
                                   /*specific_angular_momentum_ulps=*/1,
                                   /*periapsis_distance_ulps=*/0,
                                   /*apoapsis_distance_ulps=*/0);
}

}  // namespace internal_kepler_orbit
}  // namespace physics
}  // namespace principia
