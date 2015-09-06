
#include "testing_utilities/solar_system.hpp"

#include <optional.hpp>
#include <string>
#include <vector>

#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using geometry::Bivector;
using geometry::Wedge;
using physics::RelativeDegreesOfFreedom;
using physics::Body;
using quantities::SpecificAngularMomentum;
using quantities::SpecificEnergy;
using quantities::GravitationalParameter;
using quantities::Mass;
using quantities::Pow;
using quantities::Quotient;
using quantities::Speed;
using quantities::Sqrt;
using si::Radian;
using ::testing::ElementsAreArray;
using ::testing::Lt;
using ::testing::Ge;

namespace testing_utilities {

class SolarSystemTest : public testing::Test {
 protected:
  // The maximal separation of |primary| and |secondary| ignoring the influence
  // of any other bodies.
  Length SemiMajorAxis(MassiveBody const& primary_body,
                       DiscreteTrajectory<ICRFJ2000Ecliptic> const& primary,
                       MassiveBody const& secondary_body,
                       DiscreteTrajectory<ICRFJ2000Ecliptic> const& secondary) {
    GravitationalParameter const μ = primary_body.gravitational_parameter() +
                                     secondary_body.gravitational_parameter();
    RelativeDegreesOfFreedom<ICRFJ2000Ecliptic> const primary_secondary =
        primary.last().degrees_of_freedom() -
        secondary.last().degrees_of_freedom();
    Vector<Length, ICRFJ2000Ecliptic> const& r =
        primary_secondary.displacement();
    Velocity<ICRFJ2000Ecliptic> const& v = primary_secondary.velocity();
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    return -μ / (2 * ε);
  }

  // The sphere of action, or Laplace sphere (Laplace 1799) is is the region
  // where the acceleration to perturbation ratio for the secondary is larger
  // than the corresponding ratio for the primary. In this region, the 2-body
  // approximation around the secondary is better than the 2-body approximation
  // around the primary.
  Length LaplaceSphereRadiusRadius(
      MassiveBody const& primary_body,
      DiscreteTrajectory<ICRFJ2000Ecliptic> const& primary,
      MassiveBody const& secondary_body,
      DiscreteTrajectory<ICRFJ2000Ecliptic> const& secondary) {
    // Assuming secondary.mass << primary.mass.
    return SemiMajorAxis(primary_body, primary, secondary_body, secondary) *
        std::pow(secondary_body.mass() /
                 primary_body.mass(), 2.0 / 5.0);
  }

  // Tests whether |tertiary| orbits |secondary| in an orbit with excentricity
  // |excentricity| within |relative_error| and, if |primary| is not null, tests
  // that |tertiary| is within the Laplace sphere of |secondary| with respect
  // to |*primary|. If |relative_error| is greater than 1E-6, it should be tight
  // within an order of magnitude.
  void TestStronglyBoundOrbit(
      double excentricity,
      double relative_error,
      MassiveBody const& tertiary_body,
      DiscreteTrajectory<ICRFJ2000Ecliptic> const& tertiary,
      MassiveBody const& secondary_body,
      DiscreteTrajectory<ICRFJ2000Ecliptic> const& secondary,
      std::experimental::optional<MassiveBody const&> const& primary_body,
      std::experimental::optional<
          DiscreteTrajectory<ICRFJ2000Ecliptic> const&> const primary,
      std::string message) {
    GravitationalParameter const μ =
        tertiary_body.gravitational_parameter() +
        secondary_body.gravitational_parameter();
    RelativeDegreesOfFreedom<ICRFJ2000Ecliptic> const tertiary_secondary =
        tertiary.last().degrees_of_freedom() -
        secondary.last().degrees_of_freedom();
    Vector<Length, ICRFJ2000Ecliptic> const& r =
        tertiary_secondary.displacement();
    Velocity<ICRFJ2000Ecliptic> const& v = tertiary_secondary.velocity();
    Bivector<SpecificAngularMomentum, ICRFJ2000Ecliptic> const h =
        Wedge(r, v) / Radian;
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    double e = Sqrt(1 + 2 * ε * Pow<2>(h.Norm() * Radian) / Pow<2>(μ));
    EXPECT_THAT(RelativeError(excentricity, e), Lt(relative_error)) << message;
    if (relative_error > 1E-6) {
      EXPECT_THAT(RelativeError(excentricity, e),
                  Ge(relative_error / 10.0)) << message;
    }
    if (primary) {
      EXPECT_THAT(r.Norm(), Lt(LaplaceSphereRadiusRadius(*primary_body,
                                                         *primary,
                                                         secondary_body,
                                                         secondary)))
          << message;
    }
  }

  void TestStronglyBoundOrbit(
      double excentricity,
      double relative_error,
      MassiveBody const& tertiary_body,
      DiscreteTrajectory<ICRFJ2000Ecliptic> const& tertiary,
      MassiveBody const& secondary_body,
      DiscreteTrajectory<ICRFJ2000Ecliptic> const& secondary,
      std::string message) {
    TestStronglyBoundOrbit(excentricity,
                           relative_error,
                           tertiary_body,
                           tertiary,
                           secondary_body,
                           secondary,
                           std::experimental::nullopt /*tertiary_body*/,
                           std::experimental::nullopt /*tertiary*/,
                           message);
  }
};

using SolarSystemDeathTest = SolarSystemTest;

TEST_F(SolarSystemDeathTest, Parent) {
  EXPECT_DEATH({
    SolarSystem::parent(SolarSystem::kSun);
  }, "has no parent");
}

TEST_F(SolarSystemTest, Name) {
  std::vector<std::string> names;
  for (int i = SolarSystem::kSun; i <= SolarSystem::kTethys; ++i) {
    names.push_back(SolarSystem::name(i));
  }
  std::string const expected_names[] = {
      "Sun", "Jupiter", "Saturn", "Neptune", "Uranus", "Earth", "Venus", "Mars",
      "Mercury", "Ganymede", "Titan", "Callisto", "Io", "Moon", "Europa",
      "Triton", "Eris", "Pluto", "Titania", "Oberon", "Rhea", "Iapetus",
      "Charon", "Ariel", "Umbriel", "Dione", "Tethys"};
  EXPECT_THAT(names, ElementsAreArray(expected_names));
}

TEST_F(SolarSystemTest, Parent) {
  std::vector<std::string> parent_names;
  for (int i = SolarSystem::kSun + 1; i <= SolarSystem::kTethys; ++i) {
    parent_names.push_back(SolarSystem::name(SolarSystem::parent(i)));
  }
  std::string const expected_parent_names[] = {
      "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Jupiter",
      "Saturn", "Jupiter", "Jupiter", "Earth", "Jupiter", "Neptune", "Sun",
      "Sun", "Uranus", "Uranus", "Saturn", "Saturn", "Pluto", "Uranus",
      "Uranus", "Saturn", "Saturn"};
  EXPECT_THAT(parent_names, ElementsAreArray(expected_parent_names));
}

// Note(egg): We cannot call this |HierarchyAtСпутник1Launch| because gtest does
// not do a unicode-friendly stringification.  We settle for the English
// romanization.
TEST_F(SolarSystemTest, HierarchyAtSputnik1Launch) {
  auto const solar_system = SolarSystem::AtСпутник1Launch(
      SolarSystem::Accuracy::kMinorAndMajorBodies);
  auto const massive_bodies = solar_system->massive_bodies();
  auto const trajectories = solar_system->trajectories();

  auto const& sun_trajectory      = *trajectories[SolarSystem::kSun];
  auto const& jupiter_trajectory  = *trajectories[SolarSystem::kJupiter];
  auto const& saturn_trajectory   = *trajectories[SolarSystem::kSaturn];
  auto const& neptune_trajectory  = *trajectories[SolarSystem::kNeptune];
  auto const& uranus_trajectory   = *trajectories[SolarSystem::kUranus];
  auto const& earth_trajectory    = *trajectories[SolarSystem::kEarth];
  auto const& venus_trajectory    = *trajectories[SolarSystem::kVenus];
  auto const& mars_trajectory     = *trajectories[SolarSystem::kMars];
  auto const& mercury_trajectory  = *trajectories[SolarSystem::kMercury];
  auto const& ganymede_trajectory = *trajectories[SolarSystem::kGanymede];
  auto const& titan_trajectory    = *trajectories[SolarSystem::kTitan];
  auto const& callisto_trajectory = *trajectories[SolarSystem::kCallisto];
  auto const& io_trajectory       = *trajectories[SolarSystem::kIo];
  auto const& moon_trajectory     = *trajectories[SolarSystem::kMoon];
  auto const& europa_trajectory   = *trajectories[SolarSystem::kEuropa];
  auto const& triton_trajectory   = *trajectories[SolarSystem::kTriton];
  auto const& eris_trajectory     = *trajectories[SolarSystem::kEris];
  auto const& pluto_trajectory    = *trajectories[SolarSystem::kPluto];
  auto const& titania_trajectory  = *trajectories[SolarSystem::kTitania];
  auto const& oberon_trajectory   = *trajectories[SolarSystem::kOberon];
  auto const& rhea_trajectory     = *trajectories[SolarSystem::kRhea];
  auto const& iapetus_trajectory  = *trajectories[SolarSystem::kIapetus];
  auto const& charon_trajectory   = *trajectories[SolarSystem::kCharon];
  auto const& ariel_trajectory    = *trajectories[SolarSystem::kAriel];
  auto const& umbriel_trajectory  = *trajectories[SolarSystem::kUmbriel];
  auto const& dione_trajectory    = *trajectories[SolarSystem::kDione];
  auto const& tethys_trajectory   = *trajectories[SolarSystem::kTethys];

  auto const& sun      = *massive_bodies[SolarSystem::kSun];
  auto const& jupiter  = *massive_bodies[SolarSystem::kJupiter];
  auto const& saturn   = *massive_bodies[SolarSystem::kSaturn];
  auto const& neptune  = *massive_bodies[SolarSystem::kNeptune];
  auto const& uranus   = *massive_bodies[SolarSystem::kUranus];
  auto const& earth    = *massive_bodies[SolarSystem::kEarth];
  auto const& venus    = *massive_bodies[SolarSystem::kVenus];
  auto const& mars     = *massive_bodies[SolarSystem::kMars];
  auto const& mercury  = *massive_bodies[SolarSystem::kMercury];
  auto const& ganymede = *massive_bodies[SolarSystem::kGanymede];
  auto const& titan    = *massive_bodies[SolarSystem::kTitan];
  auto const& callisto = *massive_bodies[SolarSystem::kCallisto];
  auto const& io       = *massive_bodies[SolarSystem::kIo];
  auto const& moon     = *massive_bodies[SolarSystem::kMoon];
  auto const& europa   = *massive_bodies[SolarSystem::kEuropa];
  auto const& triton   = *massive_bodies[SolarSystem::kTriton];
  auto const& eris     = *massive_bodies[SolarSystem::kEris];
  auto const& pluto    = *massive_bodies[SolarSystem::kPluto];
  auto const& titania  = *massive_bodies[SolarSystem::kTitania];
  auto const& oberon   = *massive_bodies[SolarSystem::kOberon];
  auto const& rhea     = *massive_bodies[SolarSystem::kRhea];
  auto const& iapetus  = *massive_bodies[SolarSystem::kIapetus];
  auto const& charon   = *massive_bodies[SolarSystem::kCharon];
  auto const& ariel    = *massive_bodies[SolarSystem::kAriel];
  auto const& umbriel  = *massive_bodies[SolarSystem::kUmbriel];
  auto const& dione    = *massive_bodies[SolarSystem::kDione];
  auto const& tethys   = *massive_bodies[SolarSystem::kTethys];

  // Reference excentricities from HORIZONS, truncated.
  // Using centre: Sun (body centre) [500@10].
  TestStronglyBoundOrbit(4.864297E-02, 1E-6,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "jupiter");
  TestStronglyBoundOrbit(5.227008E-02, 1E-6,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "saturn");
  TestStronglyBoundOrbit(2.798871E-03, 1E-6,
                         neptune, neptune_trajectory,
                         sun, sun_trajectory,
                         "neptune");
  TestStronglyBoundOrbit(5.010917E-02, 1E-6,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "uranus");
  TestStronglyBoundOrbit(1.699349E-02, 1E-6,
                         earth, earth_trajectory,
                         sun, sun_trajectory,
                         "earth");
  TestStronglyBoundOrbit(6.797882E-03, 1E-6,
                         venus, venus_trajectory,
                         sun, sun_trajectory,
                         "venus");
  TestStronglyBoundOrbit(9.336207E-02, 1E-6,
                         mars, mars_trajectory,
                         sun, sun_trajectory,
                         "mars");
  TestStronglyBoundOrbit(2.056249E-01, 1E-6,
                         mercury, mercury_trajectory,
                         sun, sun_trajectory,
                         "mercury");
  TestStronglyBoundOrbit(2.545944E-01, 1E-6,
                         pluto, pluto_trajectory,
                         sun, sun_trajectory,
                         "pluto");
  TestStronglyBoundOrbit(4.425162E-01, 1E-6,
                         eris, eris_trajectory,
                         sun, sun_trajectory,
                         "eris");
  // Using centre: Jupiter (body centre) [500@599].
  TestStronglyBoundOrbit(2.825065E-04, 1E-4,
                         ganymede, ganymede_trajectory,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "ganymede");
  TestStronglyBoundOrbit(7.625971E-03, 1E-6,
                         callisto, callisto_trajectory,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "callisto");
  TestStronglyBoundOrbit(4.333647E-03, 1E-5,
                         io, io_trajectory,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "io");
  TestStronglyBoundOrbit(9.077806E-03, 1E-6,
                         europa, europa_trajectory,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "europa");
  // Using centre: Saturn (body centre) [500@699].
  TestStronglyBoundOrbit(2.887478E-02, 1E-6,
                         titan, titan_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "titan");
  TestStronglyBoundOrbit(8.926369E-04, 1E-5,
                         rhea, rhea_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "rhea");
  TestStronglyBoundOrbit(2.799919E-02, 1E-6,
                         iapetus, iapetus_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "iapetus");
  TestStronglyBoundOrbit(2.211120E-03, 1E-5,
                         dione, dione_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "dione");
  TestStronglyBoundOrbit(9.814475E-04, 1E-4,
                         tethys, tethys_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "tethys");
  // Using centre: Geocentric [500].
  TestStronglyBoundOrbit(5.811592E-02, 1E-6,
                         moon, moon_trajectory,
                         earth, earth_trajectory,
                         sun, sun_trajectory,
                         "moon");
  // Using centre: Neptune (body centre) [500@899]
  TestStronglyBoundOrbit(1.587851E-05, 2E-1,
                         triton, triton_trajectory,
                         neptune, neptune_trajectory,
                         sun, sun_trajectory,
                         "triton");
  // Using centre: Uranus (body centre) [500@799]
  TestStronglyBoundOrbit(1.413687E-03, 3E-3,
                         titania, titania_trajectory,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "titania");
  TestStronglyBoundOrbit(1.217327E-03, 2E-3,
                         oberon, oberon_trajectory,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "oberon");
  TestStronglyBoundOrbit(1.750702E-03, 2E-3,
                         ariel, ariel_trajectory,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "ariel");
  TestStronglyBoundOrbit(4.337777E-03, 3E-4,
                         umbriel, umbriel_trajectory,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "umbriel");
  // Using centre: Pluto (body centre) [500@999]
  TestStronglyBoundOrbit(5.077777E-05, 1E-6,
                         charon, charon_trajectory,
                         pluto, pluto_trajectory,
                         sun, sun_trajectory,
                         "charon");
}

TEST_F(SolarSystemTest, HierarchyAtSputnik2Launch) {
  auto const solar_system = SolarSystem::AtСпутник2Launch(
      SolarSystem::Accuracy::kMinorAndMajorBodies);
  auto const massive_bodies = solar_system->massive_bodies();
  auto const trajectories = solar_system->trajectories();

  auto const& sun_trajectory      = *trajectories[SolarSystem::kSun];
  auto const& jupiter_trajectory  = *trajectories[SolarSystem::kJupiter];
  auto const& saturn_trajectory   = *trajectories[SolarSystem::kSaturn];
  auto const& neptune_trajectory  = *trajectories[SolarSystem::kNeptune];
  auto const& uranus_trajectory   = *trajectories[SolarSystem::kUranus];
  auto const& earth_trajectory    = *trajectories[SolarSystem::kEarth];
  auto const& venus_trajectory    = *trajectories[SolarSystem::kVenus];
  auto const& mars_trajectory     = *trajectories[SolarSystem::kMars];
  auto const& mercury_trajectory  = *trajectories[SolarSystem::kMercury];
  auto const& ganymede_trajectory = *trajectories[SolarSystem::kGanymede];
  auto const& titan_trajectory    = *trajectories[SolarSystem::kTitan];
  auto const& callisto_trajectory = *trajectories[SolarSystem::kCallisto];
  auto const& io_trajectory       = *trajectories[SolarSystem::kIo];
  auto const& moon_trajectory     = *trajectories[SolarSystem::kMoon];
  auto const& europa_trajectory   = *trajectories[SolarSystem::kEuropa];
  auto const& triton_trajectory   = *trajectories[SolarSystem::kTriton];
  auto const& eris_trajectory     = *trajectories[SolarSystem::kEris];
  auto const& pluto_trajectory    = *trajectories[SolarSystem::kPluto];
  auto const& titania_trajectory  = *trajectories[SolarSystem::kTitania];
  auto const& oberon_trajectory   = *trajectories[SolarSystem::kOberon];
  auto const& rhea_trajectory     = *trajectories[SolarSystem::kRhea];
  auto const& iapetus_trajectory  = *trajectories[SolarSystem::kIapetus];
  auto const& charon_trajectory   = *trajectories[SolarSystem::kCharon];
  auto const& ariel_trajectory    = *trajectories[SolarSystem::kAriel];
  auto const& umbriel_trajectory  = *trajectories[SolarSystem::kUmbriel];
  auto const& dione_trajectory    = *trajectories[SolarSystem::kDione];
  auto const& tethys_trajectory   = *trajectories[SolarSystem::kTethys];

  auto const& sun      = *massive_bodies[SolarSystem::kSun];
  auto const& jupiter  = *massive_bodies[SolarSystem::kJupiter];
  auto const& saturn   = *massive_bodies[SolarSystem::kSaturn];
  auto const& neptune  = *massive_bodies[SolarSystem::kNeptune];
  auto const& uranus   = *massive_bodies[SolarSystem::kUranus];
  auto const& earth    = *massive_bodies[SolarSystem::kEarth];
  auto const& venus    = *massive_bodies[SolarSystem::kVenus];
  auto const& mars     = *massive_bodies[SolarSystem::kMars];
  auto const& mercury  = *massive_bodies[SolarSystem::kMercury];
  auto const& ganymede = *massive_bodies[SolarSystem::kGanymede];
  auto const& titan    = *massive_bodies[SolarSystem::kTitan];
  auto const& callisto = *massive_bodies[SolarSystem::kCallisto];
  auto const& io       = *massive_bodies[SolarSystem::kIo];
  auto const& moon     = *massive_bodies[SolarSystem::kMoon];
  auto const& europa   = *massive_bodies[SolarSystem::kEuropa];
  auto const& triton   = *massive_bodies[SolarSystem::kTriton];
  auto const& eris     = *massive_bodies[SolarSystem::kEris];
  auto const& pluto    = *massive_bodies[SolarSystem::kPluto];
  auto const& titania  = *massive_bodies[SolarSystem::kTitania];
  auto const& oberon   = *massive_bodies[SolarSystem::kOberon];
  auto const& rhea     = *massive_bodies[SolarSystem::kRhea];
  auto const& iapetus  = *massive_bodies[SolarSystem::kIapetus];
  auto const& charon   = *massive_bodies[SolarSystem::kCharon];
  auto const& ariel    = *massive_bodies[SolarSystem::kAriel];
  auto const& umbriel  = *massive_bodies[SolarSystem::kUmbriel];
  auto const& dione    = *massive_bodies[SolarSystem::kDione];
  auto const& tethys   = *massive_bodies[SolarSystem::kTethys];

  // Reference excentricities from HORIZONS, truncated.
  // Using centre: Sun (body centre) [500@10].
  TestStronglyBoundOrbit(4.899607E-02, 1E-6,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "jupiter");
  TestStronglyBoundOrbit(5.215911E-02, 1E-6,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "saturn");
  TestStronglyBoundOrbit(2.719093E-03, 1E-6,
                         neptune, neptune_trajectory,
                         sun, sun_trajectory,
                         "neptune");
  TestStronglyBoundOrbit(5.004209E-02, 1E-6,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "uranus");
  TestStronglyBoundOrbit(1.671840E-02, 1E-6,
                         earth, earth_trajectory,
                         sun, sun_trajectory,
                         "earth");
  TestStronglyBoundOrbit(6.792333E-03, 1E-6,
                         venus, venus_trajectory,
                         sun, sun_trajectory,
                         "venus");
  TestStronglyBoundOrbit(9.334796E-02, 1E-6,
                         mars, mars_trajectory,
                         sun, sun_trajectory,
                         "mars");
  TestStronglyBoundOrbit(2.056279E-01, 1E-6,
                         mercury, mercury_trajectory,
                         sun, sun_trajectory,
                         "mercury");
  TestStronglyBoundOrbit(2.537103E-01, 1E-6,
                         pluto, pluto_trajectory,
                         sun, sun_trajectory,
                         "pluto");
  TestStronglyBoundOrbit(4.424299E-01, 1E-6,
                         eris, eris_trajectory,
                         sun, sun_trajectory,
                         "eris");
  // Using centre: Jupiter (body centre) [500@599].
  TestStronglyBoundOrbit(4.306439E-04, 1E-4,
                         ganymede, ganymede_trajectory,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "ganymede");
  TestStronglyBoundOrbit(7.138518E-03, 1E-5,
                         callisto, callisto_trajectory,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "callisto");
  TestStronglyBoundOrbit(4.460632E-03, 1E-5,
                         io, io_trajectory,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "io");
  TestStronglyBoundOrbit(9.509972E-03, 1E-6,
                         europa, europa_trajectory,
                         jupiter, jupiter_trajectory,
                         sun, sun_trajectory,
                         "europa");
  // Using centre: Saturn (body centre) [500@699].
  TestStronglyBoundOrbit(2.882510E-02, 1E-6,
                         titan, titan_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "titan");
  TestStronglyBoundOrbit(1.228346E-03, 1E-5,
                         rhea, rhea_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "rhea");
  TestStronglyBoundOrbit(2.720904E-02, 1E-6,
                         iapetus, iapetus_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "iapetus");
  TestStronglyBoundOrbit(2.693662E-03, 1E-5,
                         dione, dione_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "dione");
  TestStronglyBoundOrbit(1.088851E-03, 1E-4,
                         tethys, tethys_trajectory,
                         saturn, saturn_trajectory,
                         sun, sun_trajectory,
                         "tethys");
  // Using centre: Geocentric [500].
  TestStronglyBoundOrbit(5.804121E-02, 1E-6,
                         moon, moon_trajectory,
                         earth, earth_trajectory,
                         sun, sun_trajectory,
                         "moon");
  // Using centre: Neptune (body centre) [500@899]
  TestStronglyBoundOrbit(1.529190E-05, 2E-1,
                         triton, triton_trajectory,
                         neptune, neptune_trajectory,
                         sun, sun_trajectory,
                         "triton");
  // Using centre: Uranus (body centre) [500@799]
  TestStronglyBoundOrbit(2.254242E-03, 3E-3,
                         titania, titania_trajectory,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "titania");
  TestStronglyBoundOrbit(4.192300E-04, 3E-3,
                         oberon, oberon_trajectory,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "oberon");
  TestStronglyBoundOrbit(2.065133E-03, 2E-3,
                         ariel, ariel_trajectory,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "ariel");
  TestStronglyBoundOrbit(3.837353E-03, 3E-4,
                         umbriel, umbriel_trajectory,
                         uranus, uranus_trajectory,
                         sun, sun_trajectory,
                         "umbriel");
  // Using centre: Pluto (body centre) [500@999]
  TestStronglyBoundOrbit(5.212037E-05, 1E-6,
                         charon, charon_trajectory,
                         pluto, pluto_trajectory,
                         sun, sun_trajectory,
                         "charon");
}

}  // namespace testing_utilities
}  // namespace principia
