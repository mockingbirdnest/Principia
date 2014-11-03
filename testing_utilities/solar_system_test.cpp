
#include "testing_utilities/solar_system.hpp"

#include <string>
#include <vector>

#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/numerics.hpp"

using principia::geometry::Bivector;
using principia::geometry::Wedge;
using principia::physics::Body;
using principia::physics::NBodySystem;
using principia::quantities::SpecificAngularMomentum;
using principia::quantities::SpecificEnergy;
using principia::quantities::GravitationalParameter;
using principia::quantities::Mass;
using principia::quantities::Pow;
using principia::quantities::Quotient;
using principia::quantities::Speed;
using principia::quantities::Sqrt;
using principia::si::Radian;
using testing::ElementsAreArray;
using testing::Lt;
using testing::Ge;

namespace principia {
namespace testing_utilities {

class SolarSystemTest : public testing::Test {
 protected:
  // The maximal separation of |primary| and |secondary| ignoring the influence
  // of any other bodies.
  Length SemiMajorAxis(Trajectory<ICRFJ2000Ecliptic> const& primary,
                       Trajectory<ICRFJ2000Ecliptic> const& secondary) {
    GravitationalParameter const μ = primary.body().gravitational_parameter() +
                                     secondary.body().gravitational_parameter();
    Vector<Length, ICRFJ2000Ecliptic> const r =
        primary.last().degrees_of_freedom().position -
        secondary.last().degrees_of_freedom().position;
    Velocity<ICRFJ2000Ecliptic> const v =
        primary.last().degrees_of_freedom().velocity -
        secondary.last().degrees_of_freedom().velocity;
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    return -μ / (2 * ε);
  }

  // The sphere of action, or Laplace sphere (Laplace 1799) is is the region
  // where the acceleration to perturbation ratio for the secondary is larger
  // than the corresponding ratio for the primary. In this region, the 2-body
  // approximation around the secondary is better than the 2-body approximation
  // around the primary.
  Length LaplaceSphereRadiusRadius(
      Trajectory<ICRFJ2000Ecliptic> const& primary,
      Trajectory<ICRFJ2000Ecliptic> const& secondary) {
    // Assuming secondary.mass << primary.mass.
    return SemiMajorAxis(primary, secondary) *
        std::pow(secondary.body().mass() / primary.body().mass(), 2.0 / 5.0);
  }

  // Tests whether |tertiary| orbits |secondary| in an orbit with excentricity
  // |excentricity| within |relative_error| and, if |primary| is not null, tests
  // that |tertiary| is within the Laplace sphere of |secondary| with respect
  // to |*primary|. If |relative_error| is greater than 1E-6, it should be tight
  // within an order of magnitude.
  void TestStronglyBoundOrbit(
      double excentricity,
      double relative_error,
      Trajectory<ICRFJ2000Ecliptic> const& tertiary,
      Trajectory<ICRFJ2000Ecliptic> const& secondary,
      Trajectory<ICRFJ2000Ecliptic> const* const primary,
      std::string message) {
    GravitationalParameter const μ = tertiary.body().gravitational_parameter() +
                                     secondary.body().gravitational_parameter();
    Vector<Length, ICRFJ2000Ecliptic> const r =
        tertiary.last().degrees_of_freedom().position -
        secondary.last().degrees_of_freedom().position;
    Velocity<ICRFJ2000Ecliptic> const v =
        tertiary.last().degrees_of_freedom().velocity -
        secondary.last().degrees_of_freedom().velocity;
    Bivector<SpecificAngularMomentum, ICRFJ2000Ecliptic> const h =
        Wedge(r, v) / Radian;
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    double e = Sqrt(1 + 2 * ε * Pow<2>(h.Norm() * Radian) / Pow<2>(μ));
    EXPECT_THAT(RelativeError(excentricity, e), Lt(relative_error)) << message;
    if (relative_error > 1E-6) {
      EXPECT_THAT(RelativeError(excentricity, e),
                  Ge(relative_error / 10.0)) << message;
    }
    if (primary != nullptr) {
      EXPECT_THAT(r.Norm(), Lt(LaplaceSphereRadiusRadius(*primary, secondary)))
          << message;
    }
  }

  std::unique_ptr<SolarSystem> solar_system_;
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
  solar_system_ = SolarSystem::AtСпутник1Launch(
      SolarSystem::Accuracy::kMinorAndMajorBodies);
  NBodySystem<ICRFJ2000Ecliptic>::Trajectories trajectories =
      solar_system_->trajectories();
  auto const& sun      = *trajectories[SolarSystem::kSun];
  auto const& jupiter  = *trajectories[SolarSystem::kJupiter];
  auto const& saturn   = *trajectories[SolarSystem::kSaturn];
  auto const& neptune  = *trajectories[SolarSystem::kNeptune];
  auto const& uranus   = *trajectories[SolarSystem::kUranus];
  auto const& earth    = *trajectories[SolarSystem::kEarth];
  auto const& venus    = *trajectories[SolarSystem::kVenus];
  auto const& mars     = *trajectories[SolarSystem::kMars];
  auto const& mercury  = *trajectories[SolarSystem::kMercury];
  auto const& ganymede = *trajectories[SolarSystem::kGanymede];
  auto const& titan    = *trajectories[SolarSystem::kTitan];
  auto const& callisto = *trajectories[SolarSystem::kCallisto];
  auto const& io       = *trajectories[SolarSystem::kIo];
  auto const& moon     = *trajectories[SolarSystem::kMoon];
  auto const& europa   = *trajectories[SolarSystem::kEuropa];
  auto const& triton   = *trajectories[SolarSystem::kTriton];
  auto const& eris     = *trajectories[SolarSystem::kEris];
  auto const& pluto    = *trajectories[SolarSystem::kPluto];
  auto const& titania  = *trajectories[SolarSystem::kTitania];
  auto const& oberon   = *trajectories[SolarSystem::kOberon];
  auto const& rhea     = *trajectories[SolarSystem::kRhea];
  auto const& iapetus  = *trajectories[SolarSystem::kIapetus];
  auto const& charon   = *trajectories[SolarSystem::kCharon];
  auto const& ariel    = *trajectories[SolarSystem::kAriel];
  auto const& umbriel  = *trajectories[SolarSystem::kUmbriel];
  auto const& dione    = *trajectories[SolarSystem::kDione];
  auto const& tethys   = *trajectories[SolarSystem::kTethys];

  // Reference excentricities from HORIZONS, truncated.
  // Using centre: Sun (body centre) [500@10].
  TestStronglyBoundOrbit(4.864297E-02, 1E-6, jupiter, sun, nullptr, "jupiter");
  TestStronglyBoundOrbit(5.227008E-02, 1E-6, saturn, sun, nullptr, "saturn");
  TestStronglyBoundOrbit(2.798871E-03, 1E-6, neptune, sun, nullptr, "neptune");
  TestStronglyBoundOrbit(5.010917E-02, 1E-6, uranus, sun, nullptr, "uranus");
  TestStronglyBoundOrbit(1.699349E-02, 1E-6, earth, sun, nullptr, "earth");
  TestStronglyBoundOrbit(6.797882E-03, 1E-6, venus, sun, nullptr, "venus");
  TestStronglyBoundOrbit(9.336207E-02, 1E-6, mars, sun, nullptr, "mars");
  TestStronglyBoundOrbit(2.056249E-01, 1E-6, mercury, sun, nullptr, "mercury");
  TestStronglyBoundOrbit(2.545944E-01, 1E-6, pluto, sun, nullptr, "pluto");
  TestStronglyBoundOrbit(4.425162E-01, 1E-6, eris, sun, nullptr, "eris");
  // Using centre: Jupiter (body centre) [500@599].
  TestStronglyBoundOrbit(2.825065E-04, 1E-4, ganymede,
                         jupiter, &sun, "ganymede");
  TestStronglyBoundOrbit(7.625971E-03, 1E-6, callisto,
                         jupiter, &sun, "callisto");
  TestStronglyBoundOrbit(4.333647E-03, 1E-5, io, jupiter, &sun, "io");
  TestStronglyBoundOrbit(9.077806E-03, 1E-6, europa, jupiter, &sun, "europa");
  // Using centre: Saturn (body centre) [500@699].
  TestStronglyBoundOrbit(2.887478E-02, 1E-6, titan, saturn, &sun, "titan");
  TestStronglyBoundOrbit(8.926369E-04, 1E-5, rhea, saturn, &sun, "rhea");
  TestStronglyBoundOrbit(2.799919E-02, 1E-6, iapetus, saturn, &sun, "iapetus");
  TestStronglyBoundOrbit(2.211120E-03, 1E-5, dione, saturn, &sun, "dione");
  TestStronglyBoundOrbit(9.814475E-04, 1E-4, tethys, saturn, &sun, "tethys");
  // Using centre: Geocentric [500].
  TestStronglyBoundOrbit(5.811592E-02, 1E-6, moon, earth, &sun, "moon");
  // Using centre: Neptune (body centre) [500@899]
  TestStronglyBoundOrbit(1.587851E-05, 2E-1, triton, neptune, &sun, "triton");
  // Using centre: Uranus (body centre) [500@799]
  TestStronglyBoundOrbit(1.413687E-03, 3E-3, titania, uranus, &sun, "titania");
  TestStronglyBoundOrbit(1.217327E-03, 2E-3, oberon, uranus, &sun, "oberon");
  TestStronglyBoundOrbit(1.750702E-03, 2E-3, ariel, uranus, &sun, "ariel");
  TestStronglyBoundOrbit(4.337777E-03, 3E-4, umbriel, uranus, &sun, "umbriel");
  // Using centre: Pluto (body centre) [500@999]
  TestStronglyBoundOrbit(5.077777E-05, 1E-6, charon, pluto, &sun, "charon");
}

TEST_F(SolarSystemTest, HierarchyAtSputnik2Launch) {
  solar_system_ = SolarSystem::AtСпутник2Launch(
      SolarSystem::Accuracy::kMinorAndMajorBodies);
  NBodySystem<ICRFJ2000Ecliptic>::Trajectories trajectories =
      solar_system_->trajectories();
  auto const& sun      = *trajectories[SolarSystem::kSun];
  auto const& jupiter  = *trajectories[SolarSystem::kJupiter];
  auto const& saturn   = *trajectories[SolarSystem::kSaturn];
  auto const& neptune  = *trajectories[SolarSystem::kNeptune];
  auto const& uranus   = *trajectories[SolarSystem::kUranus];
  auto const& earth    = *trajectories[SolarSystem::kEarth];
  auto const& venus    = *trajectories[SolarSystem::kVenus];
  auto const& mars     = *trajectories[SolarSystem::kMars];
  auto const& mercury  = *trajectories[SolarSystem::kMercury];
  auto const& ganymede = *trajectories[SolarSystem::kGanymede];
  auto const& titan    = *trajectories[SolarSystem::kTitan];
  auto const& callisto = *trajectories[SolarSystem::kCallisto];
  auto const& io       = *trajectories[SolarSystem::kIo];
  auto const& moon     = *trajectories[SolarSystem::kMoon];
  auto const& europa   = *trajectories[SolarSystem::kEuropa];
  auto const& triton   = *trajectories[SolarSystem::kTriton];
  auto const& eris     = *trajectories[SolarSystem::kEris];
  auto const& pluto    = *trajectories[SolarSystem::kPluto];
  auto const& titania  = *trajectories[SolarSystem::kTitania];
  auto const& oberon   = *trajectories[SolarSystem::kOberon];
  auto const& rhea     = *trajectories[SolarSystem::kRhea];
  auto const& iapetus  = *trajectories[SolarSystem::kIapetus];
  auto const& charon   = *trajectories[SolarSystem::kCharon];
  auto const& ariel    = *trajectories[SolarSystem::kAriel];
  auto const& umbriel  = *trajectories[SolarSystem::kUmbriel];
  auto const& dione    = *trajectories[SolarSystem::kDione];
  auto const& tethys   = *trajectories[SolarSystem::kTethys];

  // Reference excentricities from HORIZONS, truncated.
  // Using centre: Sun (body centre) [500@10].
  TestStronglyBoundOrbit(4.899607E-02, 1E-6, jupiter, sun, nullptr, "jupiter");
  TestStronglyBoundOrbit(5.215911E-02, 1E-6, saturn, sun, nullptr, "saturn");
  TestStronglyBoundOrbit(2.719093E-03, 1E-6, neptune, sun, nullptr, "neptune");
  TestStronglyBoundOrbit(5.004209E-02, 1E-6, uranus, sun, nullptr, "uranus");
  TestStronglyBoundOrbit(1.671840E-02, 1E-6, earth, sun, nullptr, "earth");
  TestStronglyBoundOrbit(6.792333E-03, 1E-6, venus, sun, nullptr, "venus");
  TestStronglyBoundOrbit(9.334796E-02, 1E-6, mars, sun, nullptr, "mars");
  TestStronglyBoundOrbit(2.056279E-01, 1E-6, mercury, sun, nullptr, "mercury");
  TestStronglyBoundOrbit(2.537103E-01, 1E-6, pluto, sun, nullptr, "pluto");
  TestStronglyBoundOrbit(4.424299E-01, 1E-6, eris, sun, nullptr, "eris");
  // Using centre: Jupiter (body centre) [500@599].
  TestStronglyBoundOrbit(4.306439E-04, 1E-4, ganymede,
                         jupiter, &sun, "ganymede");
  TestStronglyBoundOrbit(7.138518E-03, 1E-5, callisto,
                         jupiter, &sun, "callisto");
  TestStronglyBoundOrbit(4.460632E-03, 1E-5, io, jupiter, &sun, "io");
  TestStronglyBoundOrbit(9.509972E-03, 1E-6, europa, jupiter, &sun, "europa");
  // Using centre: Saturn (body centre) [500@699].
  TestStronglyBoundOrbit(2.882510E-02, 1E-6, titan, saturn, &sun, "titan");
  TestStronglyBoundOrbit(1.228346E-03, 1E-5, rhea, saturn, &sun, "rhea");
  TestStronglyBoundOrbit(2.720904E-02, 1E-6, iapetus, saturn, &sun, "iapetus");
  TestStronglyBoundOrbit(2.693662E-03, 1E-5, dione, saturn, &sun, "dione");
  TestStronglyBoundOrbit(1.088851E-03, 1E-4, tethys, saturn, &sun, "tethys");
  // Using centre: Geocentric [500].
  TestStronglyBoundOrbit(5.804121E-02, 1E-6, moon, earth, &sun, "moon");
  // Using centre: Neptune (body centre) [500@899]
  TestStronglyBoundOrbit(1.529190E-05, 2E-1, triton, neptune, &sun, "triton");
  // Using centre: Uranus (body centre) [500@799]
  TestStronglyBoundOrbit(2.254242E-03, 3E-3, titania, uranus, &sun, "titania");
  TestStronglyBoundOrbit(4.192300E-04, 3E-3, oberon, uranus, &sun, "oberon");
  TestStronglyBoundOrbit(2.065133E-03, 2E-3, ariel, uranus, &sun, "ariel");
  TestStronglyBoundOrbit(3.837353E-03, 3E-4, umbriel, uranus, &sun, "umbriel");
  // Using centre: Pluto (body centre) [500@999]
  TestStronglyBoundOrbit(5.212037E-05, 1E-6, charon, pluto, &sun, "charon");
}

}  // namespace testing_utilities
}  // namespace principia
