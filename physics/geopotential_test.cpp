#include "physics/geopotential.hpp"

#include "gtest/gtest.h"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/vanishes_before.hpp"

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
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::VanishesBefore;
using ::testing::An;

class GeopotentialTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;

  GeopotentialTest() : geopotential_(&body_) {}

  // The axis of rotation is along the z axis for ease of testing.
  AngularFrequency const angular_frequency_ = -1.5 * Radian / Second;
  Angle const right_ascension_of_pole_ = 0 * Degree;
  Angle const declination_of_pole_ = 90 * Degree;

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
  // The acceleration at a point located on the axis is along the axis.
  {
    auto const acceleration = geopotential_.SphericalHarmonicsAcceleration(
        Instant(),
        Displacement<World>({0 * Metre, 0 * Metre, 10 * Metre}),
        100 * Pow<2>(Metre),
        1.0e-3 * Pow<-3>(Metre));
    EXPECT_THAT(acceleration,
                Componentwise(VanishesBefore(1 * Pow<-2>(Metre), 0),
                              VanishesBefore(1 * Pow<-2>(Metre), 0),
                              An<Exponentiation<Length, -2>>()));
  }

  // The acceleration at a point located in the equatorial plane is directed to
  // the centre.
  {
    auto const acceleration = geopotential_.SphericalHarmonicsAcceleration(
        Instant(),
        Displacement<World>({30 * Metre, 40 * Metre, 0 * Metre}),
        2500 * Pow<2>(Metre),
        8.0e-6 * Pow<-3>(Metre));
    EXPECT_THAT(acceleration.coordinates().x / acceleration.coordinates().y,
                AlmostEquals(0.75, 1));
    EXPECT_THAT(acceleration.coordinates().z,
                VanishesBefore(1 * Pow<-2>(Metre), 0));
  }
}

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
