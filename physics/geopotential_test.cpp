#include "physics/geopotential.hpp"

#include "gtest/gtest.h"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace physics {
namespace internal_geopotential {

using geometry::Frame;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Order2ZonalCoefficient;
using quantities::Pow;
using quantities::SIUnit;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

class GeopotentialTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;

  GeopotentialTest() : geopotential_(&body_) {}

  AngularFrequency const angular_frequency_ = -1.5 * Radian / Second;
  Angle const right_ascension_of_pole_ = 37 * Degree;
  Angle const declination_of_pole_ = 123 * Degree;

  OblateBody<World> const body_ =
      OblateBody<World>(17 * SIUnit<GravitationalParameter>(),
                        RotatingBody<World>::Parameters(
                            1 * Metre,
                            3 * Radian,
                            Instant() + 4 * Second,
                            angular_frequency_,
                            right_ascension_of_pole_,
                            declination_of_pole_),
                        OblateBody<World>::Parameters(
                            163 * SIUnit<Order2ZonalCoefficient>()));

  Geopotential<World> geopotential_;
};

TEST_F(GeopotentialTest, J2) {
  auto const acceleration1 = geopotential_.SphericalHarmonicsAcceleration(
      Instant(),
      Displacement<World>({0 * Metre, 0 * Metre, 10 * Metre}),
      100 * Pow<2>(Metre),
      1.0e-3 * Pow<-3>(Metre));
}

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
