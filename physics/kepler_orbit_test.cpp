
#include "physics/kepler_orbit.hpp"

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "physics/solar_system.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace physics {
namespace internal_kepler_orbit {

using astronomy::ICRFJ2000Equator;
using astronomy::JulianDate;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;

class KeplerOrbitTest : public ::testing::Test {};

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

TEST_F(KeplerOrbitTest, EarthMoon) {
  SolarSystem<ICRFJ2000Equator> solar_system;
  solar_system.Initialize(
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
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity                = 4.772161502830355e-02;
  elements.semimajor_axis              = 3.870051955415476e+05 * Kilo(Metre);
  elements.inclination                 = 1.842335956339145e+01 * Degree;
  elements.longitude_of_ascending_node = 1.752118723367974e+00 * Degree;
  elements.argument_of_periapsis       = 3.551364385683149e+02 * Degree;
  elements.mean_anomaly                = 2.963020996150547e+02 * Degree;
  KeplerOrbit<ICRFJ2000Equator> moon_orbit(*earth, *moon, elements, date);
  Displacement<ICRFJ2000Equator> const expected_displacement(
      { 1.177367562036580e+05 * Kilo(Metre),
       -3.419908628150604e+05 * Kilo(Metre),
       -1.150659799281941e+05 * Kilo(Metre)});
  Velocity<ICRFJ2000Equator> const expected_velocity(
      {9.745048087261129e-01 * (Kilo(Metre) / Second),
       3.500672337210811e-01 * (Kilo(Metre) / Second),
       1.066306010215636e-01 * (Kilo(Metre) / Second)});
  EXPECT_THAT(moon_orbit.StateVectors(date).displacement(),
              AlmostEquals(expected_displacement, 13));
  EXPECT_THAT(moon_orbit.StateVectors(date).velocity(),
              AlmostEquals(expected_velocity, 12));
  EXPECT_THAT(*moon_orbit.elements_at_epoch().mean_motion,
              AlmostEquals(1.511718576836574e-04 * (Degree / Second), 2));

  elements.semimajor_axis = std::experimental::nullopt;
  elements.mean_motion = 1.511718576836574e-04 * (Degree / Second);
  KeplerOrbit<ICRFJ2000Equator> moon_orbit_n(*earth, *moon, elements, date);
  EXPECT_THAT(moon_orbit_n.StateVectors(date).displacement(),
              AlmostEquals(expected_displacement, 13, 15));
  EXPECT_THAT(moon_orbit_n.StateVectors(date).velocity(),
              AlmostEquals(expected_velocity, 12));

  KeplerOrbit<ICRFJ2000Equator> moon_orbit_from_state_vectors(
      *earth,
      *moon,
      {expected_displacement, expected_velocity},
      date);
  EXPECT_THAT(moon_orbit_from_state_vectors.elements_at_epoch().eccentricity,
              AlmostEquals(moon_orbit.elements_at_epoch().eccentricity, 8));
  EXPECT_THAT(*moon_orbit_from_state_vectors.elements_at_epoch().semimajor_axis,
              AlmostEquals(*moon_orbit.elements_at_epoch().semimajor_axis, 1));
  EXPECT_THAT(*moon_orbit_from_state_vectors.elements_at_epoch().mean_motion,
              AlmostEquals(*moon_orbit.elements_at_epoch().mean_motion, 1));
  EXPECT_THAT(moon_orbit_from_state_vectors.elements_at_epoch().inclination,
              AlmostEquals(moon_orbit.elements_at_epoch().inclination, 1));
  EXPECT_THAT(
      moon_orbit_from_state_vectors.elements_at_epoch()
          .longitude_of_ascending_node,
      AlmostEquals(moon_orbit.elements_at_epoch().longitude_of_ascending_node,
                   28));
  EXPECT_THAT(
      moon_orbit_from_state_vectors.elements_at_epoch().argument_of_periapsis,
      AlmostEquals(moon_orbit.elements_at_epoch().argument_of_periapsis, 6));
  EXPECT_THAT(moon_orbit_from_state_vectors.elements_at_epoch().mean_anomaly,
              AlmostEquals(moon_orbit.elements_at_epoch().mean_anomaly, 6));
}

}  // namespace internal_kepler_orbit
}  // namespace physics
}  // namespace principia
