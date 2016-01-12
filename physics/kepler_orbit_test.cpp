#include "physics/kepler_orbit.hpp"

#include "astronomy/frames.hpp"
#include "geometry/epoch.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/solar_system.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::JulianDate;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;

namespace physics {

class KeplerOrbitTest : public ::testing::Test {};

// Target body name : Moon(301) { source: DE431mx }
// Center body name : Earth(399) { source: DE431mx }
// Output units    : KM-S, deg, Julian day number (Tp)
// Reference frame : ICRF / J2000.0
// Coordinate systm : Earth Mean Equator and Equinox of Reference Epoch
// System GM : 4.0350323550225975E+05 km^3/s^2
// 2457397.500000000 = A.D. 2016-Jan-10 00:00:00.0000 (TDB)
// EC= 4.772161502830355E-02 QR= 3.685366825859605E+05 IN= 1.842335956339145E+01
// OM= 1.752118723367974E+00 W = 3.551364385683149E+02 Tp=  2457402.376861531753
// N = 1.511718576836574E-04 MA= 2.963020996150547E+02 TA= 2.912732951134421E+02
// A = 3.870051955415476E+05 AD= 4.054737084971346E+05 PR= 2.381395621619845E+06
// X = 1.177367562036580E+05 Y =-3.419908628150604E+05 Z =-1.150659799281941E+05
// VX= 9.745048087261129E-01 VY= 3.500672337210811E-01 VZ= 1.066306010215636E-01

TEST_F(KeplerOrbitTest, EarthMoon) {
  SolarSystem<ICRFJ2000Equator> solar_system;
  solar_system.Initialize(
      SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "initial_state_jd_2433282_500000000.proto.txt");
  auto const earth = SolarSystem<ICRFJ2000Equator>::MakeMassiveBody(
                         solar_system.gravity_model_message("Earth"));
  auto const moon = SolarSystem<ICRFJ2000Equator>::MakeMassiveBody(
                        solar_system.gravity_model_message("Moon"));
  // The numbers in the gravity models and those from the query above both come
  // from DE431, so the sums are the same up to round-off.
  EXPECT_THAT(
      earth->gravitational_parameter() + moon->gravitational_parameter(),
      AlmostEquals(
          4.0350323550225975E+05 * (Pow<3>(Kilo(Metre)) / Pow<2>(Second)), 1));
  Instant const date = JulianDate(2457397.500000000);
  KeplerianElements<ICRFJ2000Equator> elements;
  elements.eccentricity                = 4.772161502830355E-02;
  elements.semimajor_axis              = 3.870051955415476E+05 * Kilo(Metre);
  elements.inclination                 = 1.842335956339145E+01 * Degree;
  elements.longitude_of_ascending_node = 1.752118723367974E+00 * Degree;
  elements.argument_of_periapsis       = 3.551364385683149E+02 * Degree;
  elements.mean_anomaly                = 2.963020996150547E+02 * Degree;
  KeplerOrbit<ICRFJ2000Equator> const moon_orbit(*earth, *moon, date, elements);
  Displacement<ICRFJ2000Equator> const expected_displacement(
      { 1.177367562036580E+05 * Kilo(Metre),
       -3.419908628150604E+05 * Kilo(Metre),
       -1.150659799281941E+05 * Kilo(Metre)});
  Velocity<ICRFJ2000Equator> const expected_velocity(
      {9.745048087261129E-01 * (Kilo(Metre) / Second),
       3.500672337210811E-01 * (Kilo(Metre) / Second),
       1.066306010215636E-01 * (Kilo(Metre) / Second)});
  EXPECT_THAT(moon_orbit.PrimocentricStateVectors(date).displacement(),
              AlmostEquals(expected_displacement, 13));
  EXPECT_THAT(moon_orbit.PrimocentricStateVectors(date).velocity(),
              AlmostEquals(expected_velocity, 12));
}

TEST_F(KeplerOrbitTest, JoolSystem) {
  using KSP =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST, true>;
  MasslessBody test_particle;
  // Gravitational parameters from the KSP wiki.
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  auto const add_body = [&bodies](
      GravitationalParameter const& μ) -> MassiveBody const& {
    bodies.emplace_back(make_not_null_unique<MassiveBody>(μ));
    return *bodies.back();
  };
  auto const& sun = add_body(1.1723328E+18 * Pow<3>(Metre) / Pow<2>(Second));
  auto const& jool = add_body(2.8252800E+14 * Pow<3>(Metre) / Pow<2>(Second));
  auto const& laythe = add_body(1.9620000E+12 * Pow<3>(Metre) / Pow<2>(Second));
  auto const& vall = add_body(2.0748150E+11 * Pow<3>(Metre) / Pow<2>(Second));
  auto const& tylo = add_body(2.8252800E+12 * Pow<3>(Metre) / Pow<2>(Second));
  auto const& bop = add_body(2.4868349E+09 * Pow<3>(Metre) / Pow<2>(Second));
  auto const& pol = add_body(7.2170208E+08 * Pow<3>(Metre) / Pow<2>(Second));

  // Elements from the KSP wiki.
  KeplerianElements<KSP> jool_elements;
  jool_elements.eccentricity = 0.05;
  jool_elements.semimajor_axis = 68'773'560'320 * Metre;
  jool_elements.inclination = 1.304 * Degree;
  jool_elements.longitude_of_ascending_node = 52 * Degree;
  jool_elements.argument_of_periapsis =  0 * Degree;
  jool_elements.mean_anomaly = 0.1 * Radian;
  KeplerianElements<KSP> laythe_elements;
  laythe_elements.eccentricity = 0;
  laythe_elements.semimajor_axis = 27'184'000 * Metre;
  laythe_elements.inclination = 0 * Degree;
  laythe_elements.longitude_of_ascending_node = 0 * Degree;
  laythe_elements.argument_of_periapsis = 0 * Degree;
  laythe_elements.mean_anomaly = 3.14 * Radian;
  KeplerianElements<KSP> vall_elements;
  vall_elements.eccentricity = 0;
  vall_elements.semimajor_axis = 43'152'000 * Metre;
  vall_elements.inclination = 0 * Degree;
  vall_elements.longitude_of_ascending_node = 0 * Degree;
  vall_elements.argument_of_periapsis = 0 * Degree;
  vall_elements.mean_anomaly = 0.9 * Radian;
  KeplerianElements<KSP> tylo_elements;
  tylo_elements.eccentricity = 0;
  tylo_elements.semimajor_axis = 68'500'000 * Metre;
  tylo_elements.inclination = 0.025 * Degree;
  tylo_elements.longitude_of_ascending_node = 0 * Degree;
  tylo_elements.argument_of_periapsis = 0 * Degree;
  tylo_elements.mean_anomaly = 3.14 * Radian;
  KeplerianElements<KSP> bop_elements;
  bop_elements.eccentricity = 0.24;
  bop_elements.semimajor_axis = 128'500'000 * Metre;
  bop_elements.inclination = 15 * Degree;
  bop_elements.longitude_of_ascending_node = 10 * Degree;
  bop_elements.argument_of_periapsis = 25 * Degree;
  bop_elements.mean_anomaly = 0.9 * Radian;
  KeplerianElements<KSP> pol_elements;
  pol_elements.eccentricity = 0.17;
  pol_elements.semimajor_axis = 179'890'000 * Metre;
  pol_elements.inclination = 4.25 * Degree;
  pol_elements.longitude_of_ascending_node = 2 * Degree;
  pol_elements.argument_of_periapsis = 15 * Degree;
  pol_elements.mean_anomaly = 0.9 * Radian;

  Instant const game_epoch;
  auto const stock_jool_orbit =
      KeplerOrbit<KSP>(sun, test_particle, game_epoch, jool_elements);
  auto const stock_laythe_orbit =
      KeplerOrbit<KSP>(jool, test_particle, game_epoch, laythe_elements);
  auto const stock_vall_orbit =
      KeplerOrbit<KSP>(jool, test_particle, game_epoch, vall_elements);
  auto const stock_tylo_orbit =
      KeplerOrbit<KSP>(jool, test_particle, game_epoch, tylo_elements);
  auto const stock_bop_orbit =
      KeplerOrbit<KSP>(jool, test_particle, game_epoch, bop_elements);
  auto const stock_pol_orbit = 
      KeplerOrbit<KSP>(jool, test_particle, game_epoch, pol_elements);

  DegreesOfFreedom<KSP> const origin = {KSP::origin, Velocity<KSP>()};
  auto const stock_jool_initial_state =
      origin + stock_jool_orbit.PrimocentricStateVectors(game_epoch);
  Ephemeris<KSP> stock_ephemeris(
      std::move(bodies),
      {origin,
       stock_jool_initial_state,
       stock_jool_initial_state +
           stock_laythe_orbit.PrimocentricStateVectors(game_epoch),
       stock_jool_initial_state +
           stock_vall_orbit.PrimocentricStateVectors(game_epoch),
       stock_jool_initial_state +
           stock_tylo_orbit.PrimocentricStateVectors(game_epoch),
       stock_jool_initial_state +
           stock_bop_orbit.PrimocentricStateVectors(game_epoch),
       stock_jool_initial_state +
           stock_pol_orbit.PrimocentricStateVectors(game_epoch)},
      game_epoch,
      McLachlanAtela1992Order5Optimal<Position<KSP>>(),
      45 * Minute,
      5 * Milli(Metre));
  stock_ephemeris.Prolong(game_epoch + 500 * Day);
}

}  // namespace physics
}  // namespace principia
