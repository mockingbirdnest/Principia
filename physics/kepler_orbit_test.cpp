
#include "physics/kepler_orbit.hpp"

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
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
using astronomy::JulianDate;
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

class KeplerOrbitTest : public ::testing::Test {
 protected:
  KeplerianElements<ICRFJ2000Equator> MoonElements() const {
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

  KeplerianElements<ICRFJ2000Equator> SimpleEllipse() const {
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
        Sqrt(3) / 2 * Pow<2>(AstronomicalUnit) / (JulianYear * Radian);
    elements.periapsis_distance = 0.5 * AstronomicalUnit;
    elements.apoapsis_distance = 1.5 * AstronomicalUnit;

    elements.inclination = 1 * Degree;
    elements.longitude_of_ascending_node = 0.5 * Radian;
    elements.argument_of_periapsis = 0.5 * Radian;
    elements.longitude_of_periapsis = 1 * Radian;

    elements.true_anomaly = π / 2 * Radian;
    elements.mean_anomaly = (π / 3 - Sqrt(3) / 4) * Radian;
    elements.hyperbolic_mean_anomaly = -NaN<Angle>();
    return elements;
  }

  void ExpectConicParametersAlmostEqual(
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
      std::int64_t const apoapsis_distance_ulps) const {
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
  Instant const date = JulianDate(2457397.500000000);

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
  KeplerOrbit<ICRFJ2000Equator> moon_orbit(
      *earth, *moon, partial_elements, date);
  EXPECT_THAT(moon_orbit.StateVectors(date).displacement(),
              AlmostEquals(expected_displacement, 13));
  EXPECT_THAT(moon_orbit.StateVectors(date).velocity(),
              AlmostEquals(expected_velocity, 12));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().mean_motion,
              AlmostEquals(*MoonElements().mean_motion, 2));

  partial_elements.semimajor_axis.reset();
  partial_elements.periapsis_distance = MoonElements().periapsis_distance;
  KeplerOrbit<ICRFJ2000Equator> moon_orbit_n(
      *earth, *moon, partial_elements, date);
  EXPECT_THAT(moon_orbit_n.StateVectors(date).displacement(),
              AlmostEquals(expected_displacement, 13, 15));
  EXPECT_THAT(moon_orbit_n.StateVectors(date).velocity(),
              AlmostEquals(expected_velocity, 12));

  KeplerOrbit<ICRFJ2000Equator> moon_orbit_from_state_vectors(
      *earth,
      *moon,
      {expected_displacement, expected_velocity},
      date);
  EXPECT_THAT(*moon_orbit_from_state_vectors.elements_at_epoch().eccentricity,
              AlmostEquals(*MoonElements().eccentricity, 8));
  EXPECT_THAT(*moon_orbit_from_state_vectors.elements_at_epoch().semimajor_axis,
              AlmostEquals(*MoonElements().semimajor_axis, 1));
  EXPECT_THAT(*moon_orbit_from_state_vectors.elements_at_epoch().mean_motion,
              AlmostEquals(*MoonElements().mean_motion, 1));
  EXPECT_THAT(*moon_orbit_from_state_vectors.elements_at_epoch().period,
              AlmostEquals(*MoonElements().period, 1));
  EXPECT_THAT(
      *moon_orbit_from_state_vectors.elements_at_epoch().periapsis_distance,
      AlmostEquals(*MoonElements().periapsis_distance, 1));
  EXPECT_THAT(
      *moon_orbit_from_state_vectors.elements_at_epoch().apoapsis_distance,
      AlmostEquals(*MoonElements().apoapsis_distance, 1));
  EXPECT_THAT(moon_orbit_from_state_vectors.elements_at_epoch().inclination,
              AlmostEquals(MoonElements().inclination, 1));
  EXPECT_THAT(moon_orbit_from_state_vectors.elements_at_epoch()
                  .longitude_of_ascending_node,
              AlmostEquals(MoonElements().longitude_of_ascending_node, 28));
  EXPECT_THAT(
      *moon_orbit_from_state_vectors.elements_at_epoch().argument_of_periapsis,
      AlmostEquals(*MoonElements().argument_of_periapsis, 6));
  EXPECT_THAT(*moon_orbit_from_state_vectors.elements_at_epoch().mean_anomaly,
              AlmostEquals(*MoonElements().mean_anomaly, 6));
  EXPECT_THAT(*moon_orbit_from_state_vectors.elements_at_epoch().true_anomaly,
              AlmostEquals(*MoonElements().true_anomaly, 1));
}

TEST_F(KeplerOrbitTest, ConicFromEccentricityAndSemimajorAxis) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity = SimpleEllipse().eccentricity;
  elements.semimajor_axis = SimpleEllipse().semimajor_axis;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
 // The inputs must not change.
  EXPECT_THAT(*elements.eccentricity, Eq(*SimpleEllipse().eccentricity));
  EXPECT_THAT(*elements.semimajor_axis, Eq(*SimpleEllipse().semimajor_axis));

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

TEST_F(KeplerOrbitTest, ConicFromEccentricityAndSemiminorAxis) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity = SimpleEllipse().eccentricity;
  elements.semiminor_axis = SimpleEllipse().semiminor_axis;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.eccentricity, Eq(*SimpleEllipse().eccentricity));
  EXPECT_THAT(*elements.semiminor_axis, Eq(*SimpleEllipse().semiminor_axis));

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

TEST_F(KeplerOrbitTest, ConicFromEccentricityAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity = SimpleEllipse().eccentricity;
  elements.semilatus_rectum = SimpleEllipse().semilatus_rectum;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.eccentricity, Eq(*SimpleEllipse().eccentricity));
  EXPECT_THAT(*elements.semilatus_rectum,
              Eq(*SimpleEllipse().semilatus_rectum));

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

TEST_F(KeplerOrbitTest, ConicFromEccentricityAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity = SimpleEllipse().eccentricity;
  elements.periapsis_distance = SimpleEllipse().periapsis_distance;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.eccentricity, Eq(*SimpleEllipse().eccentricity));
  EXPECT_THAT(*elements.periapsis_distance,
              Eq(*SimpleEllipse().periapsis_distance));

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

TEST_F(KeplerOrbitTest, ConicFromEccentricityAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity = SimpleEllipse().eccentricity;
  elements.apoapsis_distance = SimpleEllipse().apoapsis_distance;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.eccentricity, Eq(*SimpleEllipse().eccentricity));
  EXPECT_THAT(*elements.apoapsis_distance,
              Eq(*SimpleEllipse().apoapsis_distance));

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

TEST_F(KeplerOrbitTest, ConicFromSemimajorAxisAndSemiminorAxis) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semimajor_axis = SimpleEllipse().semimajor_axis;
  elements.semiminor_axis = SimpleEllipse().semiminor_axis;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.semimajor_axis, Eq(*SimpleEllipse().semimajor_axis));
  EXPECT_THAT(*elements.semiminor_axis, Eq(*SimpleEllipse().semiminor_axis));

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

TEST_F(KeplerOrbitTest, ConicFromSemimajorAxisAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semimajor_axis = SimpleEllipse().semimajor_axis;
  elements.semilatus_rectum = SimpleEllipse().semilatus_rectum;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.semimajor_axis, Eq(*SimpleEllipse().semimajor_axis));
  EXPECT_THAT(*elements.semilatus_rectum,
              Eq(*SimpleEllipse().semilatus_rectum));

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

TEST_F(KeplerOrbitTest, ConicFromSemimajorAxisAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semimajor_axis = SimpleEllipse().semimajor_axis;
  elements.periapsis_distance = SimpleEllipse().periapsis_distance;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.semimajor_axis, Eq(*SimpleEllipse().semimajor_axis));
  EXPECT_THAT(*elements.periapsis_distance,
              Eq(*SimpleEllipse().periapsis_distance));

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

TEST_F(KeplerOrbitTest, ConicFromSemimajorAxisAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semimajor_axis = SimpleEllipse().semimajor_axis;
  elements.apoapsis_distance = SimpleEllipse().apoapsis_distance;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.semimajor_axis, Eq(*SimpleEllipse().semimajor_axis));
  EXPECT_THAT(*elements.apoapsis_distance,
              Eq(*SimpleEllipse().apoapsis_distance));

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

TEST_F(KeplerOrbitTest, ConicFromSemiminorAxisAndSemilatusRectum) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semiminor_axis = SimpleEllipse().semiminor_axis;
  elements.semilatus_rectum = SimpleEllipse().semilatus_rectum;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.semiminor_axis, Eq(*SimpleEllipse().semiminor_axis));
  EXPECT_THAT(*elements.semilatus_rectum,
              Eq(*SimpleEllipse().semilatus_rectum));

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

TEST_F(KeplerOrbitTest, ConicFromSemiminorAxisAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semiminor_axis = SimpleEllipse().semiminor_axis;
  elements.periapsis_distance = SimpleEllipse().periapsis_distance;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.semiminor_axis, Eq(*SimpleEllipse().semiminor_axis));
  EXPECT_THAT(*elements.periapsis_distance,
              Eq(*SimpleEllipse().periapsis_distance));

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

TEST_F(KeplerOrbitTest, ConicFromSemiminorAxisAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semiminor_axis = SimpleEllipse().semiminor_axis;
  elements.apoapsis_distance = SimpleEllipse().apoapsis_distance;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.semiminor_axis, Eq(*SimpleEllipse().semiminor_axis));
  EXPECT_THAT(*elements.apoapsis_distance,
              Eq(*SimpleEllipse().apoapsis_distance));

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

TEST_F(KeplerOrbitTest, ConicFromSemilatusRectumAndPeriapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semilatus_rectum = SimpleEllipse().semilatus_rectum;
  elements.periapsis_distance = SimpleEllipse().periapsis_distance;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.semilatus_rectum,
              Eq(*SimpleEllipse().semilatus_rectum));
  EXPECT_THAT(*elements.periapsis_distance,
              Eq(*SimpleEllipse().periapsis_distance));

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

TEST_F(KeplerOrbitTest, ConicFromSemilatusRectumAndApoapsisDistance) {
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.semilatus_rectum = SimpleEllipse().semilatus_rectum;
  elements.apoapsis_distance = SimpleEllipse().apoapsis_distance;
  // Leaving the orientation parameters and anomalies to their default values.
  // This test does not exercise them.
  elements.argument_of_periapsis.emplace();
  elements.mean_anomaly.emplace();
  elements = KeplerOrbit<ICRFJ2000Equator>(body_,
                                           MasslessBody{},
                                           elements,
                                           J2000).elements_at_epoch();
  // The inputs must not change.
  EXPECT_THAT(*elements.semilatus_rectum,
              Eq(*SimpleEllipse().semilatus_rectum));
  EXPECT_THAT(*elements.apoapsis_distance,
              Eq(*SimpleEllipse().apoapsis_distance));

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

}  // namespace internal_kepler_orbit
}  // namespace physics
}  // namespace principia
