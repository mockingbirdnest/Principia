
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
    system_ = SolarSystemAtСпутникLaunch();
  }
  // The maximal separation of |primary| and |secondary| ignoring the influence
  // of any other bodies.
  Length SemiMajorAxis(
      Body const& primary,
      Trajectory<ICRFJ2000EclipticFrame> const& primary_trajectory,
      Body const& secondary,
      Trajectory<ICRFJ2000EclipticFrame> const& secondary_trajectory) {
    GravitationalParameter const μ = primary.gravitational_parameter() +
                                     secondary.gravitational_parameter();
    Vector<Length, ICRFJ2000EclipticFrame> const r =
        primary_trajectory.Positions().rbegin()->second -
        secondary_trajectory.Positions().rbegin()->second;
    Vector<Speed, ICRFJ2000EclipticFrame> const v =
        primary_trajectory.Velocities().rbegin()->second -
        secondary_trajectory.Velocities().rbegin()->second;
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    return -μ / (2 * ε);
  }

  // The sphere of action, or Laplace sphere (Laplace 1799) is is the region
  // where the acceleration to perturbation ratio for the secondary is larger
  // than the corresponding ratio for the primary. In this region, the 2-body
  // approximation around the secondary is better than the 2-body approximation
  // around the primary.
  Length LaplaceSphereRadiusRadius(
      Body const& primary,
      Trajectory<ICRFJ2000EclipticFrame> const& primary_trajectory,
      Body const& secondary,
      Trajectory<ICRFJ2000EclipticFrame> const& secondary_trajectory) {
    // Assuming secondary.mass << primary.mass.
    return SemiMajorAxis(primary, primary_trajectory,
                         secondary, secondary_trajectory) *
        std::pow(secondary.mass() / primary.mass(), 2.0 / 5.0);
  }

  // Tests whether |tertiary| orbits |secondary| in an orbit with excentricity
  // |excentricity| within |relative_error| and, if |primary| is not null, tests
  // that |tertiary| is within the Laplace sphere of |secondary| with respect
  // to |*primary|.
  void TestStronglyBoundOrbit(
      double excentricity,
      double relative_error,
      Body const& tertiary,
      Trajectory<ICRFJ2000EclipticFrame> const& tertiary_trajectory,
      Body const& secondary,
      Trajectory<ICRFJ2000EclipticFrame> const& secondary_trajectory,
      Body const* const primary,
      Trajectory<ICRFJ2000EclipticFrame> const* const primary_trajectory,
      std::string message) {
    GravitationalParameter const μ = tertiary.gravitational_parameter() +
                                     secondary.gravitational_parameter();
    Vector<Length, ICRFJ2000EclipticFrame> const r =
        tertiary_trajectory.Positions().rbegin()->second -
        secondary_trajectory.Positions().rbegin()->second;
    Vector<Speed, ICRFJ2000EclipticFrame> const v =
        tertiary_trajectory.Velocities().rbegin()->second -
        secondary_trajectory.Velocities().rbegin()->second;
    Bivector<SpecificAngularMomentum, ICRFJ2000EclipticFrame> const h =
        Wedge(r, v) / Radian;
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    double e = Sqrt(1 + 2 * ε * Pow<2>(h.Norm() * Radian) / Pow<2>(μ));
    EXPECT_THAT(RelativeError(excentricity, e), Lt(relative_error)) << message;
    if (primary != nullptr) {
      EXPECT_THAT(r.Norm(),
                  Lt(LaplaceSphereRadiusRadius(*primary,
                                               *primary_trajectory,
                                               secondary,
                                               secondary_trajectory)))
          << message;
    }
  }

  std::unique_ptr<NBodySystem<ICRFJ2000EclipticFrame>> system_;
};

TEST_F(SolarSystemTest, Hierarchy) {
  Body const& sun      = *system_->bodies()[0];
  Body const& jupiter  = *system_->bodies()[1];
  Body const& saturn   = *system_->bodies()[2];
  Body const& neptune  = *system_->bodies()[3];
  Body const& uranus   = *system_->bodies()[4];
  Body const& earth    = *system_->bodies()[5];
  Body const& venus    = *system_->bodies()[6];
  Body const& mars     = *system_->bodies()[7];
  Body const& mercury  = *system_->bodies()[8];
  Body const& ganymede = *system_->bodies()[9];
  Body const& titan    = *system_->bodies()[10];
  Body const& callisto = *system_->bodies()[11];
  Body const& io       = *system_->bodies()[12];
  Body const& moon     = *system_->bodies()[13];
  Body const& europa   = *system_->bodies()[14];
  Body const& triton   = *system_->bodies()[15];
  Body const& eris     = *system_->bodies()[16];
  Body const& pluto    = *system_->bodies()[17];
  Trajectory<ICRFJ2000EclipticFrame> const& sun_trajectory      =
      *system_->trajectories()[0];
  Trajectory<ICRFJ2000EclipticFrame> const& jupiter_trajectory  =
      *system_->trajectories()[1];
  Trajectory<ICRFJ2000EclipticFrame> const& saturn_trajectory   =
      *system_->trajectories()[2];
  Trajectory<ICRFJ2000EclipticFrame> const& neptune_trajectory  =
      *system_->trajectories()[3];
  Trajectory<ICRFJ2000EclipticFrame> const& uranus_trajectory   =
      *system_->trajectories()[4];
  Trajectory<ICRFJ2000EclipticFrame> const& earth_trajectory    =
      *system_->trajectories()[5];
  Trajectory<ICRFJ2000EclipticFrame> const& venus_trajectory    =
      *system_->trajectories()[6];
  Trajectory<ICRFJ2000EclipticFrame> const& mars_trajectory     =
      *system_->trajectories()[7];
  Trajectory<ICRFJ2000EclipticFrame> const& mercury_trajectory  =
      *system_->trajectories()[8];
  Trajectory<ICRFJ2000EclipticFrame> const& ganymede_trajectory =
      *system_->trajectories()[9];
  Trajectory<ICRFJ2000EclipticFrame> const& titan_trajectory    =
      *system_->trajectories()[10];
  Trajectory<ICRFJ2000EclipticFrame> const& callisto_trajectory =
      *system_->trajectories()[11];
  Trajectory<ICRFJ2000EclipticFrame> const& io_trajectory       =
      *system_->trajectories()[12];
  Trajectory<ICRFJ2000EclipticFrame> const& moon_trajectory     =
      *system_->trajectories()[13];
  Trajectory<ICRFJ2000EclipticFrame> const& europa_trajectory   =
      *system_->trajectories()[14];
  Trajectory<ICRFJ2000EclipticFrame> const& triton_trajectory   =
      *system_->trajectories()[15];
  Trajectory<ICRFJ2000EclipticFrame> const& eris_trajectory     =
      *system_->trajectories()[16];
  Trajectory<ICRFJ2000EclipticFrame> const& pluto_trajectory    =
      *system_->trajectories()[17];
  // Reference excentricities from HORIZONS, truncated.
  // Using center: Sun (body center) [500@10].
  TestStronglyBoundOrbit(4.864297E-02, 1E-6, jupiter, jupiter_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "jupiter");
  TestStronglyBoundOrbit(5.227008E-02, 1E-6, saturn, saturn_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "saturn");
  TestStronglyBoundOrbit(2.798871E-03, 1E-6, neptune, neptune_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "neptune");
  TestStronglyBoundOrbit(5.010917E-02, 1E-6, uranus, uranus_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "uranus");
  TestStronglyBoundOrbit(1.699349E-02, 1E-6, earth, earth_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "earth");
  TestStronglyBoundOrbit(6.797882E-03, 1E-6, venus, venus_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "venus");
  TestStronglyBoundOrbit(9.336207E-02, 1E-6, mars, mars_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "mars");
  TestStronglyBoundOrbit(2.056249E-01, 1E-6, mercury, mercury_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "mercury");
  TestStronglyBoundOrbit(2.545944E-01, 1E-6, pluto, pluto_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "pluto");
  TestStronglyBoundOrbit(4.425162E-01, 1E-6, eris, eris_trajectory,
                         sun, sun_trajectory, nullptr, nullptr, "eris");
  // Using center: Jupiter (body center) [500@599].
  TestStronglyBoundOrbit(2.825065E-04, 1E-3, ganymede, ganymede_trajectory,
                         jupiter, jupiter_trajectory, &sun, &sun_trajectory,
                         "ganymede");
  TestStronglyBoundOrbit(7.625971E-03, 1E-4, callisto, callisto_trajectory,
                         jupiter, jupiter_trajectory, &sun, &sun_trajectory,
                         "callisto");
  TestStronglyBoundOrbit(4.333647E-03, 1E-4, io, io_trajectory, jupiter,
                         jupiter_trajectory, &sun, &sun_trajectory,
                         "io");
  TestStronglyBoundOrbit(9.077806E-03, 1E-4, europa, europa_trajectory,
                         jupiter, jupiter_trajectory, &sun, &sun_trajectory,
                         "europa");
  // Using center: Saturn (body center) [500@699].
  TestStronglyBoundOrbit(2.887478E-02, 1E-6, titan, titan_trajectory,
                         saturn, saturn_trajectory, &sun, &sun_trajectory,
                         "titan");
  // Using center: Geocentric [500].
  TestStronglyBoundOrbit(5.811592E-02, 1E-6, moon, moon_trajectory,
                         earth, earth_trajectory, &sun, &sun_trajectory,
                         "moon");
  // Using center: Neptune (body center) [500@899]
  TestStronglyBoundOrbit(1.587851E-05, 2E-1, triton, triton_trajectory,
                         neptune, neptune_trajectory, &sun, &sun_trajectory,
                         "triton");
}

}  // namespace testing_utilities
}  // namespace principia
