
#include "physics/solar_system.hpp"

#include <experimental/filesystem>

#include "astronomy/frames.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace physics {
namespace internal_solar_system {

using astronomy::ICRFJ2000Equator;
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

class SolarSystemTest : public ::testing::Test {
 protected:
  SolarSystem<ICRFJ2000Equator> solar_system_;
};

TEST_F(SolarSystemTest, RealSolarSystem) {
  solar_system_.Initialize(
      SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "initial_state_jd_2433282_500000000.proto.txt");

  EXPECT_EQ(Instant() - 50 * 365.25 * Day, solar_system_.epoch());
  EXPECT_THAT(solar_system_.names(),
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
  EXPECT_EQ(1, solar_system_.index("Callisto"));
  EXPECT_EQ(32, solar_system_.index("Venus"));

  auto const ephemeris = solar_system_.MakeEphemeris(
      /*fitting_tolerance=*/1 * Metre,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
          integrators::McLachlanAtela1992Order4Optimal<
              Position<ICRFJ2000Equator>>(),
          /*step=*/1 * Second));
  auto const earth = solar_system_.massive_body(*ephemeris, "Earth");
  EXPECT_LT(RelativeError(5.97258 * Yotta(Kilogram), earth->mass()), 6e-9);
  auto const& earth_trajectory = solar_system_.trajectory(*ephemeris, "Earth");
  EXPECT_TRUE(earth_trajectory.empty());
  EXPECT_EQ("Earth", earth->name());

  auto const& sun_initial_state = solar_system_.cartesian_initial_state_message("Sun");
  EXPECT_EQ("+1.309126697236264e+05 km", sun_initial_state.x());
  EXPECT_EQ("-7.799754996220354e-03 km/s", sun_initial_state.vx());
  auto const& sun_gravity_model = solar_system_.gravity_model_message("Sun");
  EXPECT_EQ("286.13 deg", sun_gravity_model.axis_right_ascension());
  EXPECT_EQ("63.87 deg", sun_gravity_model.axis_declination());
}

TEST_F(SolarSystemTest, KSPSystem) {
  solar_system_.Initialize(
      SOLUTION_DIR / "astronomy" / "ksp_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" / "ksp_initial_state_0_0.proto.txt");

  EXPECT_EQ(Instant{}, solar_system_.epoch());
  EXPECT_THAT(solar_system_.names(),
              ElementsAreArray({"Bop",
                                "Dres",
                                "Duna",
                                "Eeloo",
                                "Eve",
                                "Gilly",
                                "Ike",
                                "Jool",
                                "Kerbin",
                                "Kerbol",
                                "Laythe",
                                "Minmus",
                                "Moho",
                                "Mun",
                                "Pol",
                                "Tylo",
                                "Vall"}));
  EXPECT_EQ(1, solar_system_.index("Dres"));
  EXPECT_EQ(8, solar_system_.index("Kerbin"));

  auto const ephemeris = solar_system_.MakeEphemeris(
      /*fitting_tolerance=*/1 * Metre,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
          integrators::McLachlanAtela1992Order4Optimal<
              Position<ICRFJ2000Equator>>(),
          /*step=*/1 * Second));
  auto const kerbin = solar_system_.massive_body(*ephemeris, "Kerbin");
  EXPECT_LT(RelativeError(52.917061 * Zetta(Kilogram), kerbin->mass()), 6e-9);
  auto const& kerbin_trajectory =
      solar_system_.trajectory(*ephemeris, "Kerbin");
  EXPECT_TRUE(kerbin_trajectory.empty());
  EXPECT_EQ("Kerbin", kerbin->name());

  auto const& eeloo_initial_state =
      solar_system_.keplerian_initial_state_message("Eeloo");
  EXPECT_EQ(2.60000000000000009e-01,
            eeloo_initial_state.elements().eccentricity());
  EXPECT_EQ("4.00223155970064009e-08 rad /s",
            eeloo_initial_state.elements().mean_motion());
  auto const& eeloo_gravity_model =
      solar_system_.gravity_model_message("Eeloo");
  EXPECT_EQ("74410814527.049576 m^3/s^2",
            eeloo_gravity_model.gravitational_parameter());
}

TEST_F(SolarSystemTest, Clear) {
  solar_system_.Initialize(
      SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "initial_state_jd_2433282_500000000.proto.txt");
  solar_system_.RemoveMassiveBody("Io");
  solar_system_.RemoveOblateness("Sun");
  EXPECT_THAT(solar_system_.names(),
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
      solar_system_.cartesian_initial_state_message("Sun");
  EXPECT_EQ("+1.309126697236264e+05 km", sun_initial_state.x());
  EXPECT_EQ("-7.799754996220354e-03 km/s", sun_initial_state.vx());
  auto const& sun_gravity_model = solar_system_.gravity_model_message("Sun");
  EXPECT_FALSE(sun_gravity_model.has_j2());
  EXPECT_FALSE(sun_gravity_model.has_reference_radius());
}

}  // namespace internal_solar_system
}  // namespace physics
}  // namespace principia
