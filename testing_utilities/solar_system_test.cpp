
#include "testing_utilities/solar_system.hpp"

#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

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

class SolarSystemTest : public testing::Test{
protected:

  void SetUp() override {
    system_ = SolarSystemAtSputnikLaunch();
  }

  // The sphere of action, or Laplace sphere (Laplace 1799) is is the region
  // where the acceleration to perturbation ratio for the secondary is larger
  // than the corresponding ratio for the primary. In this region, the 2-body
  // approximation around the secondary is better than the 2-body approximation
  // around the primary.
  Length LaplaceSphereRadiusRadius(Body<ICRFJ2000EclipticFrame> primary,
                                   Body<ICRFJ2000EclipticFrame> secondary) {
    // Assuming a circular orbit here.
    Length const semi_major_axis = (primary.positions().back() -
                                    secondary.positions().back()).Norm();
    return semi_major_axis * std::pow(secondary.mass / primary.mass, 2.0 / 5.0);
  }

  // Tests whether |tertiary| orbits |secondary| in an orbit with excentricity
  // less than |max_excentricity| and, if |primary| is not null, tests that
  // |tertiary| is within the Laplace sphere of |secondary| with respect to |*primary|.
  void TestStronglyBoundOrbit(
      double max_excentricity,
      Body<ICRFJ2000EclipticFrame> tertiary,
      Body<ICRFJ2000EclipticFrame> secondary,
      Body<ICRFJ2000EclipticFrame>* primary = nullptr) {
    GravitationalParameter const μ = tertiary.gravitational_parameter() +
                                     secondary.gravitational_parameter();
    Vector<Length, ICRFJ2000EclipticFrame> const r =
        tertiary.positions().back() - secondary.positions().back();
    Vector<Speed, ICRFJ2000EclipticFrame> const v = 
        tertiary.velocities().back() - secondary.velocities().back();
    Bivector<SpecificAngularMomentum, ICRFJ2000EclipticFrame> const h =
        Wedge(r, v) / Radian;
    SpecificEnergy const ε = -Pow<2>(v.Norm()) / 2 - μ / r.Norm();
    double e = Sqrt(1 + 2 * ε * Pow<2>(h.Norm() * Radian) / Pow<2>(μ));
    EXPECT_THAT(e, Lt(max_excentricity));
    if (primary != nullptr) {
      EXPECT_THAT(r.Norm(), Lt(LaplaceSphereRadiusRadius(*primary, secondary)));
    }
  }

  NBodySystem<ICRFJ2000EclipticFrame> system_;
};

}
}
