
#include "physics/solar_system.hpp"

#include <ios>

#include "astronomy/frames.hpp"
#include "base/fingerprint2011.hpp"
#include "integrators/methods.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace physics {
namespace internal_solar_system {

using astronomy::ICRS;
using base::Fingerprint2011;
using geometry::Inertial;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::McLachlanAtela1992Order4Optimal;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Second;
using quantities::si::Yotta;
using quantities::si::Zetta;
using testing_utilities::RelativeError;
using ::testing::ElementsAreArray;

class SolarSystemTest : public ::testing::Test {};

TEST_F(SolarSystemTest, RealSolarSystem) {
  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2433282_500000000.proto.txt");

  EXPECT_EQ(Instant() - 50 * 365.25 * Day, solar_system.epoch());
  EXPECT_THAT(solar_system.names(),
              ElementsAreArray({"Ariel",
                                "Callisto",
                                "Ceres",
                                "Charon",
                                "Deimos",
                                "Dione",
                                "Earth",
                                "Enceladus",
                                "Eris",
                                "Europa",
                                "Ganymede",
                                "Iapetus",
                                "Io",
                                "Jupiter",
                                "Mars",
                                "Mercury",
                                "Mimas",
                                "Miranda",
                                "Moon",
                                "Neptune",
                                "Oberon",
                                "Phobos",
                                "Pluto",
                                "Rhea",
                                "Saturn",
                                "Sun",
                                "Tethys",
                                "Titan",
                                "Titania",
                                "Triton",
                                "Umbriel",
                                "Uranus",
                                "Venus",
                                "Vesta"}));
  EXPECT_EQ(1, solar_system.index("Callisto"));
  EXPECT_EQ(32, solar_system.index("Venus"));

  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Metre,
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<McLachlanAtela1992Order4Optimal,
                                                Position<ICRS>>(),
          /*step=*/1 * Second));
  auto const earth = solar_system.massive_body(*ephemeris, "Earth");
  EXPECT_LT(RelativeError(5.9723653 * Yotta(Kilogram), earth->mass()), 7e-9);
  auto const& earth_trajectory = solar_system.trajectory(*ephemeris, "Earth");
  EXPECT_TRUE(earth_trajectory.empty());
  EXPECT_EQ("Earth", earth->name());

  auto const& sun_initial_state =
      solar_system.cartesian_initial_state_message("Sun");
  EXPECT_EQ("+1.309126697236264e+05 km", sun_initial_state.x());
  EXPECT_EQ("-7.799754996220354e-03 km/s", sun_initial_state.vx());
  auto const& sun_gravity_model = solar_system.gravity_model_message("Sun");
  EXPECT_EQ("286.13 deg", sun_gravity_model.axis_right_ascension());
  EXPECT_EQ("63.87 deg", sun_gravity_model.axis_declination());
}

TEST_F(SolarSystemTest, KSPSystem) {
  using KSP = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST,
                    Inertial>;

  SolarSystem<KSP> solar_system(
      SOLUTION_DIR / "astronomy" / "kerbol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" / "kerbol_initial_state_0_0.proto.txt");

  EXPECT_EQ(Instant{}, solar_system.epoch());
  EXPECT_THAT(solar_system.names(),
              ElementsAreArray({"Bop",
                                "Dres",
                                "Duna",
                                "Eeloo",
                                "Eve",
                                "Gilly",
                                "Ike",
                                "Jool",
                                "Kerbin",
                                "Laythe",
                                "Minmus",
                                "Moho",
                                "Mun",
                                "Pol",
                                "Sun",
                                "Tylo",
                                "Vall"}));
  EXPECT_EQ(1, solar_system.index("Dres"));
  EXPECT_EQ(8, solar_system.index("Kerbin"));

  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Metre,
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<KSP>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<McLachlanAtela1992Order4Optimal,
                                                Position<KSP>>(),
          /*step=*/1 * Second));
  auto const kerbin = solar_system.massive_body(*ephemeris, "Kerbin");
  EXPECT_LT(RelativeError(52.915158 * Zetta(Kilogram), kerbin->mass()), 7e-9);
  auto const& kerbin_trajectory =
      solar_system.trajectory(*ephemeris, "Kerbin");
  EXPECT_TRUE(kerbin_trajectory.empty());
  EXPECT_EQ("Kerbin", kerbin->name());

  auto const& eeloo_initial_state =
      solar_system.keplerian_initial_state_message("Eeloo");
  EXPECT_EQ(2.60000000000000009e-01,
            eeloo_initial_state.elements().eccentricity());
  EXPECT_EQ("4.00223155970063943e-08 rad /s",
            eeloo_initial_state.elements().mean_motion());
  auto const& eeloo_gravity_model =
      solar_system.gravity_model_message("Eeloo");
  EXPECT_EQ("74410814527.0496 m^3/s^2",
            eeloo_gravity_model.gravitational_parameter());
}

TEST_F(SolarSystemTest, Clear) {
  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2433282_500000000.proto.txt");
  solar_system.RemoveMassiveBody("Io");
  solar_system.LimitOblatenessToDegree("Sun", /*max_degree=*/0);
  EXPECT_THAT(solar_system.names(),
              ElementsAreArray({"Ariel",
                                "Callisto",
                                "Ceres",
                                "Charon",
                                "Deimos",
                                "Dione",
                                "Earth",
                                "Enceladus",
                                "Eris",
                                "Europa",
                                "Ganymede",
                                "Iapetus",
                                "Jupiter",
                                "Mars",
                                "Mercury",
                                "Mimas",
                                "Miranda",
                                "Moon",
                                "Neptune",
                                "Oberon",
                                "Phobos",
                                "Pluto",
                                "Rhea",
                                "Saturn",
                                "Sun",
                                "Tethys",
                                "Titan",
                                "Titania",
                                "Triton",
                                "Umbriel",
                                "Uranus",
                                "Venus",
                                "Vesta"}));

  auto const& sun_initial_state =
      solar_system.cartesian_initial_state_message("Sun");
  EXPECT_EQ("+1.309126697236264e+05 km", sun_initial_state.x());
  EXPECT_EQ("-7.799754996220354e-03 km/s", sun_initial_state.vx());
  auto const& sun_gravity_model = solar_system.gravity_model_message("Sun");
  EXPECT_FALSE(sun_gravity_model.has_j2());
  EXPECT_FALSE(sun_gravity_model.has_reference_radius());
}

}  // namespace internal_solar_system
}  // namespace physics
}  // namespace principia
