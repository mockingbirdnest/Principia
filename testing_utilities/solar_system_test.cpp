
#include "testing_utilities/solar_system.hpp"

#include <string>

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
using testing::Lt;

namespace principia {
namespace testing_utilities {

class SolarSystemTest : public testing::Test {
 protected:
  void SetUp() {
    solar_system_ = SolarSystem::AtСпутникLaunch();
    n_body_system_.reset(new NBodySystem<ICRFJ2000Ecliptic>());
  }
  // The maximal separation of |primary| and |secondary| ignoring the influence
  // of any other bodies.
  Length SemiMajorAxis(Trajectory<ICRFJ2000Ecliptic> const& primary,
                       Trajectory<ICRFJ2000Ecliptic> const& secondary) {
    GravitationalParameter const μ = primary.body().gravitational_parameter() +
                                     secondary.body().gravitational_parameter();
    Vector<Length, ICRFJ2000Ecliptic> const r =
        primary.last_position() - secondary.last_position();
    Vector<Speed, ICRFJ2000Ecliptic> const v =
        primary.last_velocity() - secondary.last_velocity();
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
  // to |*primary|.
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
        tertiary.last_position() - secondary.last_position();
    Vector<Speed, ICRFJ2000Ecliptic> const v =
        tertiary.last_velocity() - secondary.last_velocity();
    Bivector<SpecificAngularMomentum, ICRFJ2000Ecliptic> const h =
        Wedge(r, v) / Radian;
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    double e = Sqrt(1 + 2 * ε * Pow<2>(h.Norm() * Radian) / Pow<2>(μ));
    EXPECT_THAT(RelativeError(excentricity, e), Lt(relative_error)) << message;
    if (primary != nullptr) {
      EXPECT_THAT(r.Norm(), Lt(LaplaceSphereRadiusRadius(*primary, secondary)))
          << message;
    }
  }

  std::unique_ptr<SolarSystem> solar_system_;
  std::unique_ptr<NBodySystem<ICRFJ2000Ecliptic>> n_body_system_;
};

TEST_F(SolarSystemTest, Hierarchy) {
  physics::NBodySystem<ICRFJ2000Ecliptic>::Trajectories trajectories =
      solar_system_->trajectories();
  Trajectory<ICRFJ2000Ecliptic> const& sun      = *trajectories[0];
  Trajectory<ICRFJ2000Ecliptic> const& jupiter  = *trajectories[1];
  Trajectory<ICRFJ2000Ecliptic> const& saturn   = *trajectories[2];
  Trajectory<ICRFJ2000Ecliptic> const& neptune  = *trajectories[3];
  Trajectory<ICRFJ2000Ecliptic> const& uranus   = *trajectories[4];
  Trajectory<ICRFJ2000Ecliptic> const& earth    = *trajectories[5];
  Trajectory<ICRFJ2000Ecliptic> const& venus    = *trajectories[6];
  Trajectory<ICRFJ2000Ecliptic> const& mars     = *trajectories[7];
  Trajectory<ICRFJ2000Ecliptic> const& mercury  = *trajectories[8];
  Trajectory<ICRFJ2000Ecliptic> const& ganymede = *trajectories[9];
  Trajectory<ICRFJ2000Ecliptic> const& titan    = *trajectories[10];
  Trajectory<ICRFJ2000Ecliptic> const& callisto = *trajectories[11];
  Trajectory<ICRFJ2000Ecliptic> const& io       = *trajectories[12];
  Trajectory<ICRFJ2000Ecliptic> const& moon     = *trajectories[13];
  Trajectory<ICRFJ2000Ecliptic> const& europa   = *trajectories[14];
  Trajectory<ICRFJ2000Ecliptic> const& triton   = *trajectories[15];
  Trajectory<ICRFJ2000Ecliptic> const& eris     = *trajectories[16];
  Trajectory<ICRFJ2000Ecliptic> const& pluto    = *trajectories[17];

  // Reference excentricities from HORIZONS, truncated.
  // Using center: Sun (body center) [500@10].
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
  // Using center: Jupiter (body center) [500@599].
  TestStronglyBoundOrbit(2.825065E-04, 1E-3, ganymede,
                         jupiter, &sun, "ganymede");
  TestStronglyBoundOrbit(7.625971E-03, 1E-4, callisto,
                         jupiter, &sun, "callisto");
  TestStronglyBoundOrbit(4.333647E-03, 1E-4, io, jupiter, &sun, "io");
  TestStronglyBoundOrbit(9.077806E-03, 1E-4, europa, jupiter, &sun, "europa");
  // Using center: Saturn (body center) [500@699].
  TestStronglyBoundOrbit(2.887478E-02, 1E-6, titan, saturn, &sun, "titan");
  // Using center: Geocentric [500].
  TestStronglyBoundOrbit(5.811592E-02, 1E-6, moon, earth, &sun, "moon");
  // Using center: Neptune (body center) [500@899]
  TestStronglyBoundOrbit(1.587851E-05, 2E-1, triton, neptune, &sun, "triton");
}

}  // namespace testing_utilities
}  // namespace principia
