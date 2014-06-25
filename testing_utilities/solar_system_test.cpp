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
  Length SemiMajorAxis(Body<ICRFJ2000EclipticFrame> const& primary,
                       Body<ICRFJ2000EclipticFrame> const& secondary) {
    GravitationalParameter const μ = primary.gravitational_parameter() +
                                     secondary.gravitational_parameter();
    Vector<Length, ICRFJ2000EclipticFrame> const r =
        primary.positions().back() - secondary.positions().back();
    Vector<Speed, ICRFJ2000EclipticFrame> const v =
        primary.velocities().back() - secondary.velocities().back();
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    return -μ / (2 * ε);
  }

  // The sphere of action, or Laplace sphere (Laplace 1799) is is the region
  // where the acceleration to perturbation ratio for the secondary is larger
  // than the corresponding ratio for the primary. In this region, the 2-body
  // approximation around the secondary is better than the 2-body approximation
  // around the primary.
  Length LaplaceSphereRadiusRadius(
      Body<ICRFJ2000EclipticFrame> const& primary,
      Body<ICRFJ2000EclipticFrame> const& secondary) {
    // Assuming secondary.mass << primary.mass.
    return SemiMajorAxis(primary, secondary) *
        std::pow(secondary.mass() / primary.mass(), 2.0 / 5.0);
  }

  // Tests whether |tertiary| orbits |secondary| in an orbit with excentricity
  // |excentricity| within |relative_error| and, if |primary| is not null, tests
  // that |tertiary| is within the Laplace sphere of |secondary| with respect
  // to |*primary|.
  void TestStronglyBoundOrbit(
      double excentricity,
      double relative_error,
      Body<ICRFJ2000EclipticFrame> const& tertiary,
      Body<ICRFJ2000EclipticFrame> const& secondary,
      Body<ICRFJ2000EclipticFrame> const* const primary,
      std::string message) {
    GravitationalParameter const μ = tertiary.gravitational_parameter() +
                                     secondary.gravitational_parameter();
    Vector<Length, ICRFJ2000EclipticFrame> const r =
        tertiary.positions().back() - secondary.positions().back();
    Vector<Speed, ICRFJ2000EclipticFrame> const v =
        tertiary.velocities().back() - secondary.velocities().back();
    Bivector<SpecificAngularMomentum, ICRFJ2000EclipticFrame> const h =
        Wedge(r, v) / Radian;
    SpecificEnergy const ε = Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    double e = Sqrt(1 + 2 * ε * Pow<2>(h.Norm() * Radian) / Pow<2>(μ));
    EXPECT_THAT(RelativeError(excentricity, e), Lt(relative_error)) << message;
    if (primary != nullptr) {
      EXPECT_THAT(r.Norm(),
                  Lt(LaplaceSphereRadiusRadius(*primary,
                                               secondary))) << message;
    }
  }

  std::unique_ptr<NBodySystem<ICRFJ2000EclipticFrame>> system_;
};

TEST_F(SolarSystemTest, Hierarchy) {
  Body<ICRFJ2000EclipticFrame> const& sun      = *system_->bodies()[0];
  Body<ICRFJ2000EclipticFrame> const& jupiter  = *system_->bodies()[1];
  Body<ICRFJ2000EclipticFrame> const& saturn   = *system_->bodies()[2];
  Body<ICRFJ2000EclipticFrame> const& neptune  = *system_->bodies()[3];
  Body<ICRFJ2000EclipticFrame> const& uranus   = *system_->bodies()[4];
  Body<ICRFJ2000EclipticFrame> const& earth    = *system_->bodies()[5];
  Body<ICRFJ2000EclipticFrame> const& venus    = *system_->bodies()[6];
  Body<ICRFJ2000EclipticFrame> const& mars     = *system_->bodies()[7];
  Body<ICRFJ2000EclipticFrame> const& mercury  = *system_->bodies()[8];
  Body<ICRFJ2000EclipticFrame> const& ganymede = *system_->bodies()[9];
  Body<ICRFJ2000EclipticFrame> const& titan    = *system_->bodies()[10];
  Body<ICRFJ2000EclipticFrame> const& callisto = *system_->bodies()[11];
  Body<ICRFJ2000EclipticFrame> const& io       = *system_->bodies()[12];
  Body<ICRFJ2000EclipticFrame> const& moon     = *system_->bodies()[13];
  Body<ICRFJ2000EclipticFrame> const& europa   = *system_->bodies()[14];
  Body<ICRFJ2000EclipticFrame> const& triton   = *system_->bodies()[15];
  Body<ICRFJ2000EclipticFrame> const& eris     = *system_->bodies()[16];
  Body<ICRFJ2000EclipticFrame> const& pluto    = *system_->bodies()[17];
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
  TestStronglyBoundOrbit(2.825065E-04, 1E-3, ganymede, jupiter, &sun,
                         "ganymede");
  TestStronglyBoundOrbit(7.625971E-03, 1E-4, callisto, jupiter, &sun,
                         "callisto");
  TestStronglyBoundOrbit(4.333647E-03, 1E-4, io, jupiter, &sun, "io");
  TestStronglyBoundOrbit(9.077806E-03, 1E-4, europa, jupiter, &sun, "europa");
  // Using center: Saturn (body center) [500@699].
  TestStronglyBoundOrbit(2.887478E-02, 1E-6, titan, saturn, &sun, "titan");
  // Using center: Geocentric [500].
  TestStronglyBoundOrbit(5.811592E-02, 1E-6, moon, earth, &sun, "moon");
  // Using center: Geocentric [500].
  TestStronglyBoundOrbit(5.811592E-02, 1E-6, moon, earth, &sun, "moon");
  // Using center: Neptune (body center) [500@899]
  TestStronglyBoundOrbit(1.587851E-05, 2E-1, triton, neptune, &sun, "triton");
}

}  // namespace testing_utilities
}  // namespace principia
