#include "physics/solar_system.hpp"

#include <algorithm>
#include <ios>

#include "absl/strings/str_replace.h"
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
using ::testing::ElementsAreArray;
using namespace principia::base::_fingerprint2011;
using namespace principia::geometry::_frame;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symplectic_runge_kutta_nyström_integrator;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_numerics;

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
          SymplecticRungeKuttaNyströmIntegrator<
              McLachlanAtela1992Order4Optimal,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
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
  using KSP = Frame<struct KSPTag, Inertial>;

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
          SymplecticRungeKuttaNyströmIntegrator<
              McLachlanAtela1992Order4Optimal,
              Ephemeris<KSP>::NewtonianMotionEquation>(),
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

TEST_F(SolarSystemTest, FingerprintCartesian) {
  auto gravity_model =
      ParseGravityModel(SOLUTION_DIR / "astronomy" /
                        "sol_gravity_model.proto.txt");
  auto initial_state =
      ParseInitialState(SOLUTION_DIR / "astronomy" /
                        "sol_initial_state_jd_2433282_500000000.proto.txt");

  auto const fingerprint1 =
      SolarSystem<ICRS>(gravity_model, initial_state).Fingerprint();

  // Check that the fingerprint is independent from the order of bodies.
  std::swap(*gravity_model.mutable_body(5), *gravity_model.mutable_body(7));
  std::swap(*initial_state.mutable_cartesian()->mutable_body(3),
            *initial_state.mutable_cartesian()->mutable_body(11));

  auto const fingerprint2 =
      SolarSystem<ICRS>(gravity_model, initial_state).Fingerprint();
  CHECK_EQ(fingerprint1, fingerprint2);

  // Check that the fingerprint is independent from the order of geopotential
  // rows.
  serialization::GravityModel::Body* moon_gravity_model = nullptr;
  for (auto& body : *gravity_model.mutable_body()) {
    if (body.name() == "Moon") {
      moon_gravity_model = &body;
      break;
    }
  }
  std::swap(*moon_gravity_model->mutable_geopotential()->mutable_row(2),
            *moon_gravity_model->mutable_geopotential()->mutable_row(3));

  auto const fingerprint3 =
      SolarSystem<ICRS>(gravity_model, initial_state).Fingerprint();
  CHECK_EQ(fingerprint1, fingerprint3);

  // Check that the fingerprint depends on the fields of the gravity model.
  gravity_model.mutable_body(13)->set_gravitational_parameter(
      absl::StrReplaceAll(gravity_model.body(13).gravitational_parameter(),
                          {{"km", "m"}}));

  auto const fingerprint4 =
      SolarSystem<ICRS>(gravity_model, initial_state).Fingerprint();
  CHECK_NE(fingerprint1, fingerprint4);

  // Check that the fingerprint depends on the fields of the initial state.
  initial_state.mutable_cartesian()->mutable_body(13)->set_x(
      absl::StrReplaceAll(initial_state.cartesian().body(13).x(),
                          {{"km", "m"}}));

  auto const fingerprint5 =
      SolarSystem<ICRS>(gravity_model, initial_state).Fingerprint();
  CHECK_NE(fingerprint4, fingerprint5);
}

TEST_F(SolarSystemTest, FingerprintKeplerian) {
  auto gravity_model =
      ParseGravityModel(SOLUTION_DIR / "astronomy" /
                        "kerbol_gravity_model.proto.txt");
  auto initial_state =
      ParseInitialState(SOLUTION_DIR / "astronomy" /
                        "kerbol_initial_state_0_0.proto.txt");

  auto const fingerprint1 =
      SolarSystem<ICRS>(gravity_model, initial_state).Fingerprint();

  // Check that the fingerprint is independent from the order of bodies.
  std::swap(*gravity_model.mutable_body(5), *gravity_model.mutable_body(7));
  std::swap(*initial_state.mutable_keplerian()->mutable_body(3),
            *initial_state.mutable_keplerian()->mutable_body(11));

  auto const fingerprint2 =
      SolarSystem<ICRS>(gravity_model, initial_state).Fingerprint();
  CHECK_EQ(fingerprint1, fingerprint2);

  // Check that the fingerprint depends on the fields of the gravity model.
  gravity_model.mutable_body(13)->set_gravitational_parameter(
      absl::StrReplaceAll(gravity_model.body(13).gravitational_parameter(),
                          {{"m", "km"}}));

  auto const fingerprint3 =
      SolarSystem<ICRS>(gravity_model, initial_state).Fingerprint();
  CHECK_NE(fingerprint1, fingerprint3);

  // Check that the fingerprint depends on the fields of the initial state.
  initial_state.mutable_keplerian()
      ->mutable_body(13)->mutable_elements()->set_mean_motion(
          absl::StrReplaceAll(
              initial_state.keplerian().body(13).elements().mean_motion(),
              {{"rad", "deg"}}));

  auto const fingerprint4 =
      SolarSystem<ICRS>(gravity_model, initial_state).Fingerprint();
  CHECK_NE(fingerprint3, fingerprint4);
}

}  // namespace internal_solar_system
}  // namespace physics
}  // namespace principia
