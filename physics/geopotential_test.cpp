
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
using quantities::Order3ZonalCoefficient;
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
using ::testing::Gt;
using ::testing::Lt;

class GeopotentialTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;

  GeopotentialTest()
      : massive_body_parameters_(17 * SIUnit<GravitationalParameter>()),
        rotating_body_parameters_(1 * Metre,
                                  3 * Radian,
                                  Instant() + 4 * Second,
                                  angular_frequency_,
                                  right_ascension_of_pole_,
                                  declination_of_pole_) {}

  Vector<Quotient<Acceleration, GravitationalParameter>, World>
  SphericalHarmonicsAcceleration(Geopotential<World> const& geopotential,
                                 Instant const& t,
                                 Displacement<World> const& r) {
    auto const r² = r.Norm²();
    auto const one_over_r³ = 1.0 / (r² * r.Norm());
    return geopotential.SphericalHarmonicsAcceleration(t, r, r², one_over_r³);
  }

  // The axis of rotation is along the z axis for ease of testing.
  AngularFrequency const angular_frequency_ = -1.5 * Radian / Second;
  Angle const right_ascension_of_pole_ = 0 * Degree;
  Angle const declination_of_pole_ = 90 * Degree;

  MassiveBody::Parameters const massive_body_parameters_;
  RotatingBody<World>::Parameters rotating_body_parameters_;

};

TEST_F(GeopotentialTest, J2) {
  OblateBody<World> const body =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters(
                            163 * SIUnit<Order2ZonalCoefficient>()));
  Geopotential<World> const geopotential(&body);

  // The acceleration at a point located on the axis is along the axis.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({0 * Metre, 0 * Metre, 10 * Metre}));
    EXPECT_THAT(acceleration,
                Componentwise(VanishesBefore(1 * Pow<-2>(Metre), 0),
                              VanishesBefore(1 * Pow<-2>(Metre), 0),
                              An<Exponentiation<Length, -2>>()));
  }

  // The acceleration at a point located in the equatorial plane is directed to
  // the centre.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({30 * Metre, 40 * Metre, 0 * Metre}));
    EXPECT_THAT(acceleration.coordinates().x / acceleration.coordinates().y,
                AlmostEquals(0.75, 1));
    EXPECT_THAT(acceleration.coordinates().z,
                VanishesBefore(1 * Pow<-2>(Metre), 0));
  }

  // The acceleration at a random point nudges the overall force away from the
  // centre and towards the equatorial plane.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({1e2 * Metre, 0 * Metre, 1e2 * Metre}));
    EXPECT_THAT(acceleration.coordinates().x, Gt(0 * Pow<-2>(Metre)));
    EXPECT_THAT(acceleration.coordinates().z, Lt(0 * Pow<-2>(Metre)));
  }
}

TEST_F(GeopotentialTest, J3) {
  OblateBody<World> const body =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters(
                            163 * SIUnit<Order2ZonalCoefficient>(),
                            -51 * SIUnit<Order3ZonalCoefficient>()));
  Geopotential<World> const geopotential(&body);

  // The acceleration at a point located on the axis is along the axis.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({0 * Metre, 0 * Metre, 10 * Metre}));
    EXPECT_THAT(acceleration,
                Componentwise(VanishesBefore(1 * Pow<-2>(Metre), 0),
                              VanishesBefore(1 * Pow<-2>(Metre), 0),
                              An<Exponentiation<Length, -2>>()));
  }

  // The acceleration at a point located in the equatorial plane points towards
  // the north, as it does on Earth (I think).
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({30 * Metre, 40 * Metre, 0 * Metre}));
    EXPECT_THAT(acceleration.coordinates().x / acceleration.coordinates().y,
                AlmostEquals(0.75, 1));
    EXPECT_THAT(acceleration.coordinates().z,
                Not(VanishesBefore(1 * Pow<-2>(Metre), 0)));
    EXPECT_THAT(acceleration.coordinates().z, Gt(0 * Pow<-2>(Metre)));
  }
}

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
