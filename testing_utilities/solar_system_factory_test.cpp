
#include "testing_utilities/solar_system_factory.hpp"

#include "base/macros.hpp"
#include OPTIONAL_HEADER
#include <string>
#include <vector>

#include "astronomy/frames.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/massive_body.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::Bivector;
using geometry::Wedge;
using physics::Body;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using physics::RelativeDegreesOfFreedom;
using quantities::SpecificAngularMomentum;
using quantities::SpecificEnergy;
using quantities::GravitationalParameter;
using quantities::Mass;
using quantities::Pow;
using quantities::Quotient;
using quantities::Speed;
using quantities::Sqrt;
using quantities::si::Radian;
using ::testing::ElementsAreArray;
using ::testing::Lt;
using ::testing::Ge;

namespace testing_utilities {

class SolarSystemFactoryTest : public testing::Test {
 protected:
  // The maximal separation of |primary| and |secondary| ignoring the influence
  // of any other bodies.
  Length SemiMajorAxis(
      MassiveBody const& primary_body,
      DegreesOfFreedom<ICRFJ2000Equator> const& primary_dof,
      MassiveBody const& secondary_body,
      DegreesOfFreedom<ICRFJ2000Equator> const& secondary_dof) {
    GravitationalParameter const μ = primary_body.gravitational_parameter() +
                                     secondary_body.gravitational_parameter();
    RelativeDegreesOfFreedom<ICRFJ2000Equator> const primary_secondary =
        primary_dof - secondary_dof;
    Vector<Length, ICRFJ2000Equator> const& r =
        primary_secondary.displacement();
    Velocity<ICRFJ2000Equator> const& v = primary_secondary.velocity();
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
      DegreesOfFreedom<ICRFJ2000Equator> const& primary_dof,
      MassiveBody const& secondary_body,
      DegreesOfFreedom<ICRFJ2000Equator> const& secondary_dof) {
    // Assuming secondary.mass << primary.mass.
    return SemiMajorAxis(primary_body, primary_dof,
                         secondary_body, secondary_dof) *
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
      DegreesOfFreedom<ICRFJ2000Equator> const& tertiary_dof,
      MassiveBody const& secondary_body,
      DegreesOfFreedom<ICRFJ2000Equator> const& secondary_dof,
      std::experimental::optional<MassiveBody const&> const& primary_body,
      std::experimental::optional<
          DegreesOfFreedom<ICRFJ2000Equator> const&> const primary_dof,
      std::string message) {
    GravitationalParameter const μ =
        tertiary_body.gravitational_parameter() +
        secondary_body.gravitational_parameter();
    RelativeDegreesOfFreedom<ICRFJ2000Equator> const tertiary_secondary =
        tertiary_dof - secondary_dof;
    Vector<Length, ICRFJ2000Equator> const& r =
        tertiary_secondary.displacement();
    Velocity<ICRFJ2000Equator> const& v = tertiary_secondary.velocity();
    Bivector<SpecificAngularMomentum, ICRFJ2000Equator> const h =
        Wedge(r, v) / Radian;
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    double e = Sqrt(1 + 2 * ε * Pow<2>(h.Norm() * Radian) / Pow<2>(μ));
    EXPECT_THAT(RelativeError(excentricity, e), Lt(relative_error)) << message;
    if (relative_error > 1E-6) {
      EXPECT_THAT(RelativeError(excentricity, e),
                  Ge(relative_error / 10.0)) << message;
    }
    if (primary_dof) {
      EXPECT_THAT(r.Norm(), Lt(LaplaceSphereRadiusRadius(*primary_body,
                                                         *primary_dof,
                                                         secondary_body,
                                                         secondary_dof)))
          << message;
    }
  }

  void TestStronglyBoundOrbit(
      double excentricity,
      double relative_error,
      MassiveBody const& tertiary_body,
      DegreesOfFreedom<ICRFJ2000Equator> const& tertiary_dof,
      MassiveBody const& secondary_body,
      DegreesOfFreedom<ICRFJ2000Equator> const& secondary_dof,
      std::string message) {
    TestStronglyBoundOrbit(excentricity,
                           relative_error,
                           tertiary_body,
                           tertiary_dof,
                           secondary_body,
                           secondary_dof,
                           std::experimental::nullopt /*tertiary_body*/,
                           std::experimental::nullopt /*tertiary*/,
                           message);
  }

  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> GetDegreesOfFreedom(
      SolarSystem<ICRFJ2000Equator> const& solar_system) {
    std::vector<DegreesOfFreedom<ICRFJ2000Equator>> degrees_of_freedom;
    for (int i = SolarSystemFactory::kSun;
         i <= SolarSystemFactory::kLastBody;
         ++i) {
      degrees_of_freedom.emplace_back(
          solar_system.initial_state(SolarSystemFactory::name(i)));
    }
    return degrees_of_freedom;
  }

  std::vector<std::unique_ptr<MassiveBody>> GetMassiveBodies(
    SolarSystem<ICRFJ2000Equator> const& solar_system) {
    std::vector<std::unique_ptr<MassiveBody>> massive_bodies;
    for (int i = SolarSystemFactory::kSun;
         i <= SolarSystemFactory::kLastBody;
         ++i) {
      massive_bodies.emplace_back(
          SolarSystem<ICRFJ2000Equator>::MakeMassiveBody(
              solar_system.gravity_model_message(SolarSystemFactory::name(i))));
    }
    return massive_bodies;
  }
};

using SolarSystemFactoryDeathTest = SolarSystemFactoryTest;

TEST_F(SolarSystemFactoryDeathTest, Parent) {
  EXPECT_DEATH({
    SolarSystemFactory::parent(SolarSystemFactory::kSun);
  }, "has no parent");
}

TEST_F(SolarSystemFactoryTest, Name) {
  std::vector<std::string> names;
  for (int i = SolarSystemFactory::kSun;
       i <= SolarSystemFactory::kLastBody;
       ++i) {
    names.push_back(SolarSystemFactory::name(i));
  }
  std::string const expected_names[] = {
      "Sun", "Jupiter", "Saturn", "Neptune", "Uranus", "Earth", "Venus", "Mars",
      "Mercury", "Ganymede", "Titan", "Callisto", "Io", "Moon", "Europa",
      "Triton", "Eris", "Pluto", "Titania", "Oberon", "Rhea", "Iapetus",
      "Charon", "Ariel", "Umbriel", "Dione", "Tethys"};
  EXPECT_THAT(names, ElementsAreArray(expected_names));
}

TEST_F(SolarSystemFactoryTest, Parent) {
  std::vector<std::string> parent_names;
  for (int i = SolarSystemFactory::kSun + 1;
       i <= SolarSystemFactory::kLastBody;
       ++i) {
    parent_names.push_back(
        SolarSystemFactory::name(SolarSystemFactory::parent(i)));
  }
  std::string const expected_parent_names[] = {
      "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Jupiter",
      "Saturn", "Jupiter", "Jupiter", "Earth", "Jupiter", "Neptune", "Sun",
      "Sun", "Uranus", "Uranus", "Saturn", "Saturn", "Pluto", "Uranus",
      "Uranus", "Saturn", "Saturn"};
  EXPECT_THAT(parent_names, ElementsAreArray(expected_parent_names));
}

TEST_F(SolarSystemFactoryTest, HierarchyAtСпутник1Launch) {
  auto const solar_system = SolarSystemFactory::AtСпутник1Launch(
      SolarSystemFactory::Accuracy::kMinorAndMajorBodies);
  auto const massive_bodies = GetMassiveBodies(*solar_system);
  auto const dof = GetDegreesOfFreedom(*solar_system);

  auto const& sun_dof      = dof[SolarSystemFactory::kSun];
  auto const& jupiter_dof  = dof[SolarSystemFactory::kJupiter];
  auto const& saturn_dof   = dof[SolarSystemFactory::kSaturn];
  auto const& neptune_dof  = dof[SolarSystemFactory::kNeptune];
  auto const& uranus_dof   = dof[SolarSystemFactory::kUranus];
  auto const& earth_dof    = dof[SolarSystemFactory::kEarth];
  auto const& venus_dof    = dof[SolarSystemFactory::kVenus];
  auto const& mars_dof     = dof[SolarSystemFactory::kMars];
  auto const& mercury_dof  = dof[SolarSystemFactory::kMercury];
  auto const& ganymede_dof = dof[SolarSystemFactory::kGanymede];
  auto const& titan_dof    = dof[SolarSystemFactory::kTitan];
  auto const& callisto_dof = dof[SolarSystemFactory::kCallisto];
  auto const& io_dof       = dof[SolarSystemFactory::kIo];
  auto const& moon_dof     = dof[SolarSystemFactory::kMoon];
  auto const& europa_dof   = dof[SolarSystemFactory::kEuropa];
  auto const& triton_dof   = dof[SolarSystemFactory::kTriton];
  auto const& eris_dof     = dof[SolarSystemFactory::kEris];
  auto const& pluto_dof    = dof[SolarSystemFactory::kPluto];
  auto const& titania_dof  = dof[SolarSystemFactory::kTitania];
  auto const& oberon_dof   = dof[SolarSystemFactory::kOberon];
  auto const& rhea_dof     = dof[SolarSystemFactory::kRhea];
  auto const& iapetus_dof  = dof[SolarSystemFactory::kIapetus];
  auto const& charon_dof   = dof[SolarSystemFactory::kCharon];
  auto const& ariel_dof    = dof[SolarSystemFactory::kAriel];
  auto const& umbriel_dof  = dof[SolarSystemFactory::kUmbriel];
  auto const& dione_dof    = dof[SolarSystemFactory::kDione];
  auto const& tethys_dof   = dof[SolarSystemFactory::kTethys];

  auto const& sun      = *massive_bodies[SolarSystemFactory::kSun];
  auto const& jupiter  = *massive_bodies[SolarSystemFactory::kJupiter];
  auto const& saturn   = *massive_bodies[SolarSystemFactory::kSaturn];
  auto const& neptune  = *massive_bodies[SolarSystemFactory::kNeptune];
  auto const& uranus   = *massive_bodies[SolarSystemFactory::kUranus];
  auto const& earth    = *massive_bodies[SolarSystemFactory::kEarth];
  auto const& venus    = *massive_bodies[SolarSystemFactory::kVenus];
  auto const& mars     = *massive_bodies[SolarSystemFactory::kMars];
  auto const& mercury  = *massive_bodies[SolarSystemFactory::kMercury];
  auto const& ganymede = *massive_bodies[SolarSystemFactory::kGanymede];
  auto const& titan    = *massive_bodies[SolarSystemFactory::kTitan];
  auto const& callisto = *massive_bodies[SolarSystemFactory::kCallisto];
  auto const& io       = *massive_bodies[SolarSystemFactory::kIo];
  auto const& moon     = *massive_bodies[SolarSystemFactory::kMoon];
  auto const& europa   = *massive_bodies[SolarSystemFactory::kEuropa];
  auto const& triton   = *massive_bodies[SolarSystemFactory::kTriton];
  auto const& eris     = *massive_bodies[SolarSystemFactory::kEris];
  auto const& pluto    = *massive_bodies[SolarSystemFactory::kPluto];
  auto const& titania  = *massive_bodies[SolarSystemFactory::kTitania];
  auto const& oberon   = *massive_bodies[SolarSystemFactory::kOberon];
  auto const& rhea     = *massive_bodies[SolarSystemFactory::kRhea];
  auto const& iapetus  = *massive_bodies[SolarSystemFactory::kIapetus];
  auto const& charon   = *massive_bodies[SolarSystemFactory::kCharon];
  auto const& ariel    = *massive_bodies[SolarSystemFactory::kAriel];
  auto const& umbriel  = *massive_bodies[SolarSystemFactory::kUmbriel];
  auto const& dione    = *massive_bodies[SolarSystemFactory::kDione];
  auto const& tethys   = *massive_bodies[SolarSystemFactory::kTethys];

  // Reference excentricities from HORIZONS, truncated.
  // Using centre: Sun (body centre) [500@10].
  TestStronglyBoundOrbit(4.864297E-02, 1E-6,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "jupiter");
  TestStronglyBoundOrbit(5.227008E-02, 1E-6,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "saturn");
  TestStronglyBoundOrbit(2.798871E-03, 1E-6,
                         neptune, neptune_dof,
                         sun, sun_dof,
                         "neptune");
  TestStronglyBoundOrbit(5.010917E-02, 1E-6,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "uranus");
  TestStronglyBoundOrbit(1.699349E-02, 1E-6,
                         earth, earth_dof,
                         sun, sun_dof,
                         "earth");
  TestStronglyBoundOrbit(6.797882E-03, 1E-6,
                         venus, venus_dof,
                         sun, sun_dof,
                         "venus");
  TestStronglyBoundOrbit(9.336207E-02, 1E-6,
                         mars, mars_dof,
                         sun, sun_dof,
                         "mars");
  TestStronglyBoundOrbit(2.056249E-01, 1E-6,
                         mercury, mercury_dof,
                         sun, sun_dof,
                         "mercury");
  TestStronglyBoundOrbit(2.545944E-01, 1E-6,
                         pluto, pluto_dof,
                         sun, sun_dof,
                         "pluto");
  TestStronglyBoundOrbit(4.425162E-01, 1E-5,
                         eris, eris_dof,
                         sun, sun_dof,
                         "eris");
  // Using centre: Jupiter (body centre) [500@599].
  TestStronglyBoundOrbit(2.825065E-04, 1E-5,
                         ganymede, ganymede_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "ganymede");
  TestStronglyBoundOrbit(7.625971E-03, 1E-6,
                         callisto, callisto_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "callisto");
  TestStronglyBoundOrbit(4.333647E-03, 1E-5,
                         io, io_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "io");
  TestStronglyBoundOrbit(9.077806E-03, 1E-6,
                         europa, europa_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "europa");
  // Using centre: Saturn (body centre) [500@699].
  TestStronglyBoundOrbit(2.887478E-02, 1E-6,
                         titan, titan_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "titan");
  TestStronglyBoundOrbit(8.926369E-04, 1E-4,
                         rhea, rhea_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "rhea");
  TestStronglyBoundOrbit(2.799919E-02, 1E-5,
                         iapetus, iapetus_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "iapetus");
  TestStronglyBoundOrbit(2.211120E-03, 1E-3,
                         dione, dione_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "dione");
  TestStronglyBoundOrbit(9.814475E-04, 1E-3,
                         tethys, tethys_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "tethys");
  // Using centre: Geocentric [500].
  TestStronglyBoundOrbit(5.811592E-02, 1E-6,
                         moon, moon_dof,
                         earth, earth_dof,
                         sun, sun_dof,
                         "moon");
  // Using centre: Neptune (body centre) [500@899]
  TestStronglyBoundOrbit(1.587851E-05, 1E-5,
                         triton, triton_dof,
                         neptune, neptune_dof,
                         sun, sun_dof,
                         "triton");
  // Using centre: Uranus (body centre) [500@799]
  TestStronglyBoundOrbit(1.413687E-03, 1E-6,
                         titania, titania_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "titania");
  TestStronglyBoundOrbit(1.217327E-03, 1E-6,
                         oberon, oberon_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "oberon");
  TestStronglyBoundOrbit(1.750702E-03, 1E-6,
                         ariel, ariel_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "ariel");
  TestStronglyBoundOrbit(4.337777E-03, 1E-6,
                         umbriel, umbriel_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "umbriel");
  // Using centre: Pluto (body centre) [500@999]
  TestStronglyBoundOrbit(5.077777E-05, 1E-6,
                         charon, charon_dof,
                         pluto, pluto_dof,
                         sun, sun_dof,
                         "charon");
}

TEST_F(SolarSystemFactoryTest, HierarchyAtСпутник2Launch) {
  auto const solar_system = SolarSystemFactory::AtСпутник2Launch(
      SolarSystemFactory::Accuracy::kMinorAndMajorBodies);
  auto const massive_bodies = GetMassiveBodies(*solar_system);
  auto const dof = GetDegreesOfFreedom(*solar_system);

  auto const& sun_dof      = dof[SolarSystemFactory::kSun];
  auto const& jupiter_dof  = dof[SolarSystemFactory::kJupiter];
  auto const& saturn_dof   = dof[SolarSystemFactory::kSaturn];
  auto const& neptune_dof  = dof[SolarSystemFactory::kNeptune];
  auto const& uranus_dof   = dof[SolarSystemFactory::kUranus];
  auto const& earth_dof    = dof[SolarSystemFactory::kEarth];
  auto const& venus_dof    = dof[SolarSystemFactory::kVenus];
  auto const& mars_dof     = dof[SolarSystemFactory::kMars];
  auto const& mercury_dof  = dof[SolarSystemFactory::kMercury];
  auto const& ganymede_dof = dof[SolarSystemFactory::kGanymede];
  auto const& titan_dof    = dof[SolarSystemFactory::kTitan];
  auto const& callisto_dof = dof[SolarSystemFactory::kCallisto];
  auto const& io_dof       = dof[SolarSystemFactory::kIo];
  auto const& moon_dof     = dof[SolarSystemFactory::kMoon];
  auto const& europa_dof   = dof[SolarSystemFactory::kEuropa];
  auto const& triton_dof   = dof[SolarSystemFactory::kTriton];
  auto const& eris_dof     = dof[SolarSystemFactory::kEris];
  auto const& pluto_dof    = dof[SolarSystemFactory::kPluto];
  auto const& titania_dof  = dof[SolarSystemFactory::kTitania];
  auto const& oberon_dof   = dof[SolarSystemFactory::kOberon];
  auto const& rhea_dof     = dof[SolarSystemFactory::kRhea];
  auto const& iapetus_dof  = dof[SolarSystemFactory::kIapetus];
  auto const& charon_dof   = dof[SolarSystemFactory::kCharon];
  auto const& ariel_dof    = dof[SolarSystemFactory::kAriel];
  auto const& umbriel_dof  = dof[SolarSystemFactory::kUmbriel];
  auto const& dione_dof    = dof[SolarSystemFactory::kDione];
  auto const& tethys_dof   = dof[SolarSystemFactory::kTethys];

  auto const& sun      = *massive_bodies[SolarSystemFactory::kSun];
  auto const& jupiter  = *massive_bodies[SolarSystemFactory::kJupiter];
  auto const& saturn   = *massive_bodies[SolarSystemFactory::kSaturn];
  auto const& neptune  = *massive_bodies[SolarSystemFactory::kNeptune];
  auto const& uranus   = *massive_bodies[SolarSystemFactory::kUranus];
  auto const& earth    = *massive_bodies[SolarSystemFactory::kEarth];
  auto const& venus    = *massive_bodies[SolarSystemFactory::kVenus];
  auto const& mars     = *massive_bodies[SolarSystemFactory::kMars];
  auto const& mercury  = *massive_bodies[SolarSystemFactory::kMercury];
  auto const& ganymede = *massive_bodies[SolarSystemFactory::kGanymede];
  auto const& titan    = *massive_bodies[SolarSystemFactory::kTitan];
  auto const& callisto = *massive_bodies[SolarSystemFactory::kCallisto];
  auto const& io       = *massive_bodies[SolarSystemFactory::kIo];
  auto const& moon     = *massive_bodies[SolarSystemFactory::kMoon];
  auto const& europa   = *massive_bodies[SolarSystemFactory::kEuropa];
  auto const& triton   = *massive_bodies[SolarSystemFactory::kTriton];
  auto const& eris     = *massive_bodies[SolarSystemFactory::kEris];
  auto const& pluto    = *massive_bodies[SolarSystemFactory::kPluto];
  auto const& titania  = *massive_bodies[SolarSystemFactory::kTitania];
  auto const& oberon   = *massive_bodies[SolarSystemFactory::kOberon];
  auto const& rhea     = *massive_bodies[SolarSystemFactory::kRhea];
  auto const& iapetus  = *massive_bodies[SolarSystemFactory::kIapetus];
  auto const& charon   = *massive_bodies[SolarSystemFactory::kCharon];
  auto const& ariel    = *massive_bodies[SolarSystemFactory::kAriel];
  auto const& umbriel  = *massive_bodies[SolarSystemFactory::kUmbriel];
  auto const& dione    = *massive_bodies[SolarSystemFactory::kDione];
  auto const& tethys   = *massive_bodies[SolarSystemFactory::kTethys];

  // Reference excentricities from HORIZONS, truncated.
  // Using centre: Sun (body centre) [500@10].
  TestStronglyBoundOrbit(4.899607E-02, 1E-6,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "jupiter");
  TestStronglyBoundOrbit(5.215911E-02, 1E-6,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "saturn");
  TestStronglyBoundOrbit(2.719093E-03, 1E-6,
                         neptune, neptune_dof,
                         sun, sun_dof,
                         "neptune");
  TestStronglyBoundOrbit(5.004209E-02, 1E-6,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "uranus");
  TestStronglyBoundOrbit(1.671840E-02, 1E-6,
                         earth, earth_dof,
                         sun, sun_dof,
                         "earth");
  TestStronglyBoundOrbit(6.792333E-03, 1E-6,
                         venus, venus_dof,
                         sun, sun_dof,
                         "venus");
  TestStronglyBoundOrbit(9.334796E-02, 1E-6,
                         mars, mars_dof,
                         sun, sun_dof,
                         "mars");
  TestStronglyBoundOrbit(2.056279E-01, 1E-6,
                         mercury, mercury_dof,
                         sun, sun_dof,
                         "mercury");
  TestStronglyBoundOrbit(2.537103E-01, 1E-6,
                         pluto, pluto_dof,
                         sun, sun_dof,
                         "pluto");
  TestStronglyBoundOrbit(4.424299E-01, 1E-5,
                         eris, eris_dof,
                         sun, sun_dof,
                         "eris");
  // Using centre: Jupiter (body centre) [500@599].
  TestStronglyBoundOrbit(4.306439E-04, 1E-5,
                         ganymede, ganymede_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "ganymede");
  TestStronglyBoundOrbit(7.138518E-03, 1E-6,
                         callisto, callisto_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "callisto");
  TestStronglyBoundOrbit(4.460632E-03, 1E-5,
                         io, io_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "io");
  TestStronglyBoundOrbit(9.509972E-03, 1E-6,
                         europa, europa_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "europa");
  // Using centre: Saturn (body centre) [500@699].
  TestStronglyBoundOrbit(2.882510E-02, 1E-6,
                         titan, titan_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "titan");
  TestStronglyBoundOrbit(1.228346E-03, 1E-4,
                         rhea, rhea_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "rhea");
  TestStronglyBoundOrbit(2.720904E-02, 1E-5,
                         iapetus, iapetus_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "iapetus");
  TestStronglyBoundOrbit(2.693662E-03, 1E-3,
                         dione, dione_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "dione");
  TestStronglyBoundOrbit(1.088851E-03, 1E-3,
                         tethys, tethys_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "tethys");
  // Using centre: Geocentric [500].
  TestStronglyBoundOrbit(5.804121E-02, 1E-6,
                         moon, moon_dof,
                         earth, earth_dof,
                         sun, sun_dof,
                         "moon");
  // Using centre: Neptune (body centre) [500@899]
  TestStronglyBoundOrbit(1.529190E-05, 1E-5,
                         triton, triton_dof,
                         neptune, neptune_dof,
                         sun, sun_dof,
                         "triton");
  // Using centre: Uranus (body centre) [500@799]
  TestStronglyBoundOrbit(2.254242E-03, 1E-6,
                         titania, titania_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "titania");
  TestStronglyBoundOrbit(4.192300E-04, 1E-6,
                         oberon, oberon_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "oberon");
  TestStronglyBoundOrbit(2.065133E-03, 1E-6,
                         ariel, ariel_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "ariel");
  TestStronglyBoundOrbit(3.837353E-03, 1E-6,
                         umbriel, umbriel_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "umbriel");
  // Using centre: Pluto (body centre) [500@999]
  TestStronglyBoundOrbit(5.212037E-05, 1E-6,
                         charon, charon_dof,
                         pluto, pluto_dof,
                         sun, sun_dof,
                         "charon");
}

}  // namespace testing_utilities
}  // namespace principia
