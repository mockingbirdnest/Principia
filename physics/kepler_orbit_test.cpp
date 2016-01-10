#include "physics/kepler_orbit.hpp"

#include "astronomy/frames.hpp"
#include "geometry/epoch.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/solar_system.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::JulianDate;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using testing_utilities::AbsoluteError;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;

namespace physics {

class KeplerOrbitTest : public ::testing::Test {};

// Target body name : Moon(301) { source: DE431mx }
// Center body name : Earth(399) { source: DE431mx }
// Output units : AU - D, deg, Julian day number(Tp)
// Reference frame : ICRF / J2000.0
// Coordinate systm : Earth Mean Equator and Equinox of Reference Epoch
// System GM : 3.8897282712664292E+05 km^3/s^2
// 2457397.500000000 = A.D. 2016 - Jan - 10 00 : 00 : 00.0000 (TDB)
// EC= 4.772161502830354E-02 QR= 3.640587465676356E+05 IN= 1.842335956339145E+01
// OM= 1.752118723367971E+00 W = 3.551364385683146E+02 Tp= 2457402.376861531753
// N = 1.511718576836573E-04 MA= 2.963020996150549E+02 TA= 2.912732951134424E+02
// A = 3.823028563002153E+05 AD= 4.005469660327950E+05 PR= 2.381395621619846E+06
// X = 1.177367562036580E+05 Y =-3.419908628150604E+05 Z =-1.150659799281941E+05
// VX= 9.745048087261129E-01 VY= 3.500672337210811E-01 VZ= 1.066306010215636E-01

TEST_F(KeplerOrbitTest, EarthMoon) {
  SolarSystem<ICRFJ2000Equator> solar_system;
  solar_system.Initialize(
      SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "initial_state_jd_2433282_500000000.proto.txt");
  auto const earth = solar_system.MakeMassiveBody(
                         solar_system.gravity_model_message("Earth"));
  auto const moon = solar_system.MakeMassiveBody(
                        solar_system.gravity_model_message("Moon"));
  EXPECT_EQ(earth->gravitational_parameter() + moon->gravitational_parameter(),
            3.8897282712664292E+05 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  TwoBodySystem const earth_moon_system{earth.get(), moon.get()};
  Instant const date = JulianDate(2457397.500000000);
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity                = 4.772161502830354E-02;
  elements.semimajor_axis              = 3.823028563002153E+05 * Kilo(Metre);
  elements.inclination                 = 1.842335956339145E+01 * Degree;
  elements.longitude_of_ascending_node = 1.752118723367971E+00 * Degree;
  elements.argument_of_periapsis       = 3.551364385683146E+02 * Degree;
  elements.mean_anomaly                = 2.963020996150549E+02 * Degree;
  KeplerOrbit<ICRFJ2000Equator> const moon_orbit(earth_moon_system,
                                                 date,
                                                 elements);
  LOG(ERROR) << moon_orbit.StateVectors(date);
}

}  // namespace physics
}  // namespace principia
