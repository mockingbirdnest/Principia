
#include "testing_utilities/solar_system_factory.hpp"

#include <optional>
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
using astronomy::J2000;
using geometry::Bivector;
using geometry::Vector;
using geometry::Velocity;
using geometry::Wedge;
using physics::Body;
using physics::DegreesOfFreedom;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::RelativeDegreesOfFreedom;
using physics::SolarSystem;
using quantities::SpecificAngularMomentum;
using quantities::SpecificEnergy;
using quantities::GravitationalParameter;
using quantities::Length;
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
    Length const a = *KeplerOrbit<ICRFJ2000Equator>{
        primary_body,
        secondary_body,
        secondary_dof - primary_dof, J2000}.elements_at_epoch().semimajor_axis;
    // Assuming secondary.mass << primary.mass.
    return a * std::pow(secondary_body.mass() / primary_body.mass(), 2.0 / 5.0);
  }

  // Tests whether |tertiary| orbits |secondary| in an orbit with excentricity
  // |excentricity| within |relative_error| and, if |primary| is not null, tests
  // that |tertiary| is within the Laplace sphere of |secondary| with respect
  // to |*primary|. If |relative_error| is greater than 1e-6, it should be tight
  // within an order of magnitude.
  void TestStronglyBoundOrbit(
      double eccentricity,
      double relative_error,
      MassiveBody const& tertiary_body,
      DegreesOfFreedom<ICRFJ2000Equator> const& tertiary_dof,
      MassiveBody const& secondary_body,
      DegreesOfFreedom<ICRFJ2000Equator> const& secondary_dof,
      std::optional<std::reference_wrapper<MassiveBody const>> const&
          primary_body,
      std::optional <std::reference_wrapper<
          DegreesOfFreedom<ICRFJ2000Equator> const>> const& primary_dof,
      std::string const& message) {
    RelativeDegreesOfFreedom<ICRFJ2000Equator> const tertiary_secondary =
        tertiary_dof - secondary_dof;
    KeplerOrbit<ICRFJ2000Equator> orbit{
        secondary_body, tertiary_body, tertiary_secondary, J2000};
    Vector<Length, ICRFJ2000Equator> const& r =
        tertiary_secondary.displacement();
    EXPECT_THAT(
        RelativeError(eccentricity, *orbit.elements_at_epoch().eccentricity),
        Lt(relative_error))
        << message;
    if (relative_error > 1e-6) {
      EXPECT_THAT(
          RelativeError(eccentricity, *orbit.elements_at_epoch().eccentricity),
          Ge(relative_error / 10.0))
          << message;
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
      std::string const& message) {
    TestStronglyBoundOrbit(excentricity,
                           relative_error,
                           tertiary_body,
                           tertiary_dof,
                           secondary_body,
                           secondary_dof,
                           /*primary_body=*/std::nullopt,
                           /*primary_dof=*/std::nullopt,
                           message);
  }

  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> GetDegreesOfFreedom(
      SolarSystem<ICRFJ2000Equator> const& solar_system) {
    std::vector<DegreesOfFreedom<ICRFJ2000Equator>> degrees_of_freedom;
    for (int i = SolarSystemFactory::Sun;
         i <= SolarSystemFactory::LastBody;
         ++i) {
      degrees_of_freedom.emplace_back(
          solar_system.degrees_of_freedom(SolarSystemFactory::name(i)));
    }
    return degrees_of_freedom;
  }

  std::vector<std::unique_ptr<MassiveBody>> GetMassiveBodies(
    SolarSystem<ICRFJ2000Equator> const& solar_system) {
    std::vector<std::unique_ptr<MassiveBody>> massive_bodies;
    for (int i = SolarSystemFactory::Sun;
         i <= SolarSystemFactory::LastBody;
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
    SolarSystemFactory::parent(SolarSystemFactory::Sun);
  }, "has no parent");
}

TEST_F(SolarSystemFactoryTest, Name) {
  std::vector<std::string> names;
  for (int i = SolarSystemFactory::Sun;
       i <= SolarSystemFactory::LastBody;
       ++i) {
    names.push_back(SolarSystemFactory::name(i));
  }
  std::string const expected_names[] = {
      "Sun", "Jupiter", "Saturn", "Neptune", "Uranus", "Earth", "Venus", "Mars",
      "Mercury", "Ganymede", "Titan", "Callisto", "Io", "Moon", "Europa",
      "Triton", "Eris", "Pluto", "Titania", "Oberon", "Rhea", "Iapetus",
      "Charon", "Ariel", "Umbriel", "Dione", "Ceres", "Tethys", "Vesta",
      "Enceladus", "Miranda", "Mimas", "Phobos", "Deimos"};
  EXPECT_THAT(names, ElementsAreArray(expected_names));
}

TEST_F(SolarSystemFactoryTest, Parent) {
  std::vector<std::string> parent_names;
  for (int i = SolarSystemFactory::Sun + 1;
       i <= SolarSystemFactory::LastBody;
       ++i) {
    parent_names.push_back(
        SolarSystemFactory::name(SolarSystemFactory::parent(i)));
  }
  std::string const expected_parent_names[] = {
      "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Jupiter",
      "Saturn", "Jupiter", "Jupiter", "Earth", "Jupiter", "Neptune", "Sun",
      "Sun", "Uranus", "Uranus", "Saturn", "Saturn", "Pluto", "Uranus",
      "Uranus", "Saturn", "Sun", "Saturn", "Sun", "Saturn", "Uranus", "Saturn",
      "Mars", "Mars"};
  EXPECT_THAT(parent_names, ElementsAreArray(expected_parent_names));
}

TEST_F(SolarSystemFactoryTest, HierarchyAtСпутник1Launch) {
  auto const solar_system = SolarSystemFactory::AtСпутник1Launch(
      SolarSystemFactory::Accuracy::MinorAndMajorBodies);
  auto const massive_bodies = GetMassiveBodies(*solar_system);
  auto const dof = GetDegreesOfFreedom(*solar_system);

  auto const& sun_dof      = dof[SolarSystemFactory::Sun];
  auto const& jupiter_dof  = dof[SolarSystemFactory::Jupiter];
  auto const& saturn_dof   = dof[SolarSystemFactory::Saturn];
  auto const& neptune_dof  = dof[SolarSystemFactory::Neptune];
  auto const& uranus_dof   = dof[SolarSystemFactory::Uranus];
  auto const& earth_dof    = dof[SolarSystemFactory::Earth];
  auto const& venus_dof    = dof[SolarSystemFactory::Venus];
  auto const& mars_dof     = dof[SolarSystemFactory::Mars];
  auto const& mercury_dof  = dof[SolarSystemFactory::Mercury];
  auto const& ganymede_dof = dof[SolarSystemFactory::Ganymede];
  auto const& titan_dof    = dof[SolarSystemFactory::Titan];
  auto const& callisto_dof = dof[SolarSystemFactory::Callisto];
  auto const& io_dof       = dof[SolarSystemFactory::Io];
  auto const& moon_dof     = dof[SolarSystemFactory::Moon];
  auto const& europa_dof   = dof[SolarSystemFactory::Europa];
  auto const& triton_dof   = dof[SolarSystemFactory::Triton];
  auto const& eris_dof     = dof[SolarSystemFactory::Eris];
  auto const& pluto_dof    = dof[SolarSystemFactory::Pluto];
  auto const& titania_dof  = dof[SolarSystemFactory::Titania];
  auto const& oberon_dof   = dof[SolarSystemFactory::Oberon];
  auto const& rhea_dof     = dof[SolarSystemFactory::Rhea];
  auto const& iapetus_dof  = dof[SolarSystemFactory::Iapetus];
  auto const& charon_dof   = dof[SolarSystemFactory::Charon];
  auto const& ariel_dof    = dof[SolarSystemFactory::Ariel];
  auto const& umbriel_dof  = dof[SolarSystemFactory::Umbriel];
  auto const& dione_dof    = dof[SolarSystemFactory::Dione];
  auto const& tethys_dof   = dof[SolarSystemFactory::Tethys];

  auto const& sun      = *massive_bodies[SolarSystemFactory::Sun];
  auto const& jupiter  = *massive_bodies[SolarSystemFactory::Jupiter];
  auto const& saturn   = *massive_bodies[SolarSystemFactory::Saturn];
  auto const& neptune  = *massive_bodies[SolarSystemFactory::Neptune];
  auto const& uranus   = *massive_bodies[SolarSystemFactory::Uranus];
  auto const& earth    = *massive_bodies[SolarSystemFactory::Earth];
  auto const& venus    = *massive_bodies[SolarSystemFactory::Venus];
  auto const& mars     = *massive_bodies[SolarSystemFactory::Mars];
  auto const& mercury  = *massive_bodies[SolarSystemFactory::Mercury];
  auto const& ganymede = *massive_bodies[SolarSystemFactory::Ganymede];
  auto const& titan    = *massive_bodies[SolarSystemFactory::Titan];
  auto const& callisto = *massive_bodies[SolarSystemFactory::Callisto];
  auto const& io       = *massive_bodies[SolarSystemFactory::Io];
  auto const& moon     = *massive_bodies[SolarSystemFactory::Moon];
  auto const& europa   = *massive_bodies[SolarSystemFactory::Europa];
  auto const& triton   = *massive_bodies[SolarSystemFactory::Triton];
  auto const& eris     = *massive_bodies[SolarSystemFactory::Eris];
  auto const& pluto    = *massive_bodies[SolarSystemFactory::Pluto];
  auto const& titania  = *massive_bodies[SolarSystemFactory::Titania];
  auto const& oberon   = *massive_bodies[SolarSystemFactory::Oberon];
  auto const& rhea     = *massive_bodies[SolarSystemFactory::Rhea];
  auto const& iapetus  = *massive_bodies[SolarSystemFactory::Iapetus];
  auto const& charon   = *massive_bodies[SolarSystemFactory::Charon];
  auto const& ariel    = *massive_bodies[SolarSystemFactory::Ariel];
  auto const& umbriel  = *massive_bodies[SolarSystemFactory::Umbriel];
  auto const& dione    = *massive_bodies[SolarSystemFactory::Dione];
  auto const& tethys   = *massive_bodies[SolarSystemFactory::Tethys];

  // Reference excentricities from HORIZONS, truncated.
  // Using centre: Sun (body centre) [500@10].
  TestStronglyBoundOrbit(4.864297e-02, 1e-6,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "jupiter");
  TestStronglyBoundOrbit(5.227008e-02, 1e-6,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "saturn");
  TestStronglyBoundOrbit(2.798871e-03, 1e-6,
                         neptune, neptune_dof,
                         sun, sun_dof,
                         "neptune");
  TestStronglyBoundOrbit(5.010917e-02, 1e-6,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "uranus");
  TestStronglyBoundOrbit(1.699349e-02, 1e-6,
                         earth, earth_dof,
                         sun, sun_dof,
                         "earth");
  TestStronglyBoundOrbit(6.797882e-03, 1e-6,
                         venus, venus_dof,
                         sun, sun_dof,
                         "venus");
  TestStronglyBoundOrbit(9.336207e-02, 1e-6,
                         mars, mars_dof,
                         sun, sun_dof,
                         "mars");
  TestStronglyBoundOrbit(2.056249e-01, 1e-6,
                         mercury, mercury_dof,
                         sun, sun_dof,
                         "mercury");
  TestStronglyBoundOrbit(2.545944e-01, 1e-6,
                         pluto, pluto_dof,
                         sun, sun_dof,
                         "pluto");
  TestStronglyBoundOrbit(4.425162e-01, 1e-5,
                         eris, eris_dof,
                         sun, sun_dof,
                         "eris");
  // Using centre: Jupiter (body centre) [500@599].
  TestStronglyBoundOrbit(2.825065e-04, 1e-5,
                         ganymede, ganymede_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "ganymede");
  TestStronglyBoundOrbit(7.625971e-03, 1e-6,
                         callisto, callisto_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "callisto");
  TestStronglyBoundOrbit(4.333647e-03, 1e-5,
                         io, io_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "io");
  TestStronglyBoundOrbit(9.077806e-03, 1e-6,
                         europa, europa_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "europa");
  // Using centre: Saturn (body centre) [500@699].
  TestStronglyBoundOrbit(2.887478e-02, 1e-6,
                         titan, titan_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "titan");
  TestStronglyBoundOrbit(8.926369e-04, 1e-4,
                         rhea, rhea_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "rhea");
  TestStronglyBoundOrbit(2.799919e-02, 1e-5,
                         iapetus, iapetus_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "iapetus");
  TestStronglyBoundOrbit(2.211120e-03, 1e-3,
                         dione, dione_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "dione");
  TestStronglyBoundOrbit(9.814475e-04, 1e-3,
                         tethys, tethys_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "tethys");
  // Using centre: Geocentric [500].
  TestStronglyBoundOrbit(5.811592e-02, 1e-6,
                         moon, moon_dof,
                         earth, earth_dof,
                         sun, sun_dof,
                         "moon");
  // Using centre: Neptune (body centre) [500@899]
  TestStronglyBoundOrbit(1.587851e-05, 1e-5,
                         triton, triton_dof,
                         neptune, neptune_dof,
                         sun, sun_dof,
                         "triton");
  // Using centre: Uranus (body centre) [500@799]
  TestStronglyBoundOrbit(1.413687e-03, 1e-6,
                         titania, titania_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "titania");
  TestStronglyBoundOrbit(1.217327e-03, 1e-6,
                         oberon, oberon_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "oberon");
  TestStronglyBoundOrbit(1.750702e-03, 1e-6,
                         ariel, ariel_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "ariel");
  TestStronglyBoundOrbit(4.337777e-03, 1e-6,
                         umbriel, umbriel_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "umbriel");
  // Using centre: Pluto (body centre) [500@999]
  TestStronglyBoundOrbit(5.077777e-05, 1e-6,
                         charon, charon_dof,
                         pluto, pluto_dof,
                         sun, sun_dof,
                         "charon");
}

TEST_F(SolarSystemFactoryTest, HierarchyAtСпутник2Launch) {
  auto const solar_system = SolarSystemFactory::AtСпутник2Launch(
      SolarSystemFactory::Accuracy::MinorAndMajorBodies);
  auto const massive_bodies = GetMassiveBodies(*solar_system);
  auto const dof = GetDegreesOfFreedom(*solar_system);

  auto const& sun_dof      = dof[SolarSystemFactory::Sun];
  auto const& jupiter_dof  = dof[SolarSystemFactory::Jupiter];
  auto const& saturn_dof   = dof[SolarSystemFactory::Saturn];
  auto const& neptune_dof  = dof[SolarSystemFactory::Neptune];
  auto const& uranus_dof   = dof[SolarSystemFactory::Uranus];
  auto const& earth_dof    = dof[SolarSystemFactory::Earth];
  auto const& venus_dof    = dof[SolarSystemFactory::Venus];
  auto const& mars_dof     = dof[SolarSystemFactory::Mars];
  auto const& mercury_dof  = dof[SolarSystemFactory::Mercury];
  auto const& ganymede_dof = dof[SolarSystemFactory::Ganymede];
  auto const& titan_dof    = dof[SolarSystemFactory::Titan];
  auto const& callisto_dof = dof[SolarSystemFactory::Callisto];
  auto const& io_dof       = dof[SolarSystemFactory::Io];
  auto const& moon_dof     = dof[SolarSystemFactory::Moon];
  auto const& europa_dof   = dof[SolarSystemFactory::Europa];
  auto const& triton_dof   = dof[SolarSystemFactory::Triton];
  auto const& eris_dof     = dof[SolarSystemFactory::Eris];
  auto const& pluto_dof    = dof[SolarSystemFactory::Pluto];
  auto const& titania_dof  = dof[SolarSystemFactory::Titania];
  auto const& oberon_dof   = dof[SolarSystemFactory::Oberon];
  auto const& rhea_dof     = dof[SolarSystemFactory::Rhea];
  auto const& iapetus_dof  = dof[SolarSystemFactory::Iapetus];
  auto const& charon_dof   = dof[SolarSystemFactory::Charon];
  auto const& ariel_dof    = dof[SolarSystemFactory::Ariel];
  auto const& umbriel_dof  = dof[SolarSystemFactory::Umbriel];
  auto const& dione_dof    = dof[SolarSystemFactory::Dione];
  auto const& tethys_dof   = dof[SolarSystemFactory::Tethys];

  auto const& sun      = *massive_bodies[SolarSystemFactory::Sun];
  auto const& jupiter  = *massive_bodies[SolarSystemFactory::Jupiter];
  auto const& saturn   = *massive_bodies[SolarSystemFactory::Saturn];
  auto const& neptune  = *massive_bodies[SolarSystemFactory::Neptune];
  auto const& uranus   = *massive_bodies[SolarSystemFactory::Uranus];
  auto const& earth    = *massive_bodies[SolarSystemFactory::Earth];
  auto const& venus    = *massive_bodies[SolarSystemFactory::Venus];
  auto const& mars     = *massive_bodies[SolarSystemFactory::Mars];
  auto const& mercury  = *massive_bodies[SolarSystemFactory::Mercury];
  auto const& ganymede = *massive_bodies[SolarSystemFactory::Ganymede];
  auto const& titan    = *massive_bodies[SolarSystemFactory::Titan];
  auto const& callisto = *massive_bodies[SolarSystemFactory::Callisto];
  auto const& io       = *massive_bodies[SolarSystemFactory::Io];
  auto const& moon     = *massive_bodies[SolarSystemFactory::Moon];
  auto const& europa   = *massive_bodies[SolarSystemFactory::Europa];
  auto const& triton   = *massive_bodies[SolarSystemFactory::Triton];
  auto const& eris     = *massive_bodies[SolarSystemFactory::Eris];
  auto const& pluto    = *massive_bodies[SolarSystemFactory::Pluto];
  auto const& titania  = *massive_bodies[SolarSystemFactory::Titania];
  auto const& oberon   = *massive_bodies[SolarSystemFactory::Oberon];
  auto const& rhea     = *massive_bodies[SolarSystemFactory::Rhea];
  auto const& iapetus  = *massive_bodies[SolarSystemFactory::Iapetus];
  auto const& charon   = *massive_bodies[SolarSystemFactory::Charon];
  auto const& ariel    = *massive_bodies[SolarSystemFactory::Ariel];
  auto const& umbriel  = *massive_bodies[SolarSystemFactory::Umbriel];
  auto const& dione    = *massive_bodies[SolarSystemFactory::Dione];
  auto const& tethys   = *massive_bodies[SolarSystemFactory::Tethys];

  // Reference excentricities from HORIZONS, truncated.
  // Using centre: Sun (body centre) [500@10].
  TestStronglyBoundOrbit(4.899607e-02, 1e-6,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "jupiter");
  TestStronglyBoundOrbit(5.215911e-02, 1e-6,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "saturn");
  TestStronglyBoundOrbit(2.719093e-03, 1e-6,
                         neptune, neptune_dof,
                         sun, sun_dof,
                         "neptune");
  TestStronglyBoundOrbit(5.004209e-02, 1e-6,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "uranus");
  TestStronglyBoundOrbit(1.671840e-02, 1e-6,
                         earth, earth_dof,
                         sun, sun_dof,
                         "earth");
  TestStronglyBoundOrbit(6.792333e-03, 1e-6,
                         venus, venus_dof,
                         sun, sun_dof,
                         "venus");
  TestStronglyBoundOrbit(9.334796e-02, 1e-6,
                         mars, mars_dof,
                         sun, sun_dof,
                         "mars");
  TestStronglyBoundOrbit(2.056279e-01, 1e-6,
                         mercury, mercury_dof,
                         sun, sun_dof,
                         "mercury");
  TestStronglyBoundOrbit(2.537103e-01, 1e-6,
                         pluto, pluto_dof,
                         sun, sun_dof,
                         "pluto");
  TestStronglyBoundOrbit(4.424299e-01, 1e-5,
                         eris, eris_dof,
                         sun, sun_dof,
                         "eris");
  // Using centre: Jupiter (body centre) [500@599].
  TestStronglyBoundOrbit(4.306439e-04, 1e-5,
                         ganymede, ganymede_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "ganymede");
  TestStronglyBoundOrbit(7.138518e-03, 1e-6,
                         callisto, callisto_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "callisto");
  TestStronglyBoundOrbit(4.460632e-03, 1e-5,
                         io, io_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "io");
  TestStronglyBoundOrbit(9.509972e-03, 1e-6,
                         europa, europa_dof,
                         jupiter, jupiter_dof,
                         sun, sun_dof,
                         "europa");
  // Using centre: Saturn (body centre) [500@699].
  TestStronglyBoundOrbit(2.882510e-02, 1e-6,
                         titan, titan_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "titan");
  TestStronglyBoundOrbit(1.228346e-03, 1e-4,
                         rhea, rhea_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "rhea");
  TestStronglyBoundOrbit(2.720904e-02, 1e-5,
                         iapetus, iapetus_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "iapetus");
  TestStronglyBoundOrbit(2.693662e-03, 1e-3,
                         dione, dione_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "dione");
  TestStronglyBoundOrbit(1.088851e-03, 1e-3,
                         tethys, tethys_dof,
                         saturn, saturn_dof,
                         sun, sun_dof,
                         "tethys");
  // Using centre: Geocentric [500].
  TestStronglyBoundOrbit(5.804121e-02, 1e-6,
                         moon, moon_dof,
                         earth, earth_dof,
                         sun, sun_dof,
                         "moon");
  // Using centre: Neptune (body centre) [500@899]
  TestStronglyBoundOrbit(1.529190e-05, 1e-6,
                         triton, triton_dof,
                         neptune, neptune_dof,
                         sun, sun_dof,
                         "triton");
  // Using centre: Uranus (body centre) [500@799]
  TestStronglyBoundOrbit(2.254242e-03, 1e-6,
                         titania, titania_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "titania");
  TestStronglyBoundOrbit(4.192300e-04, 1e-6,
                         oberon, oberon_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "oberon");
  TestStronglyBoundOrbit(2.065133e-03, 1e-6,
                         ariel, ariel_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "ariel");
  TestStronglyBoundOrbit(3.837353e-03, 1e-6,
                         umbriel, umbriel_dof,
                         uranus, uranus_dof,
                         sun, sun_dof,
                         "umbriel");
  // Using centre: Pluto (body centre) [500@999]
  TestStronglyBoundOrbit(5.212037e-05, 1e-6,
                         charon, charon_dof,
                         pluto, pluto_dof,
                         sun, sun_dof,
                         "charon");
}

}  // namespace testing_utilities
}  // namespace principia
