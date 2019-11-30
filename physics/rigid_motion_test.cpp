
#include "physics/rigid_motion.hpp"

#include "geometry/frame.hpp"
#include "geometry/permutation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {
namespace internal_rigid_motion {

using geometry::AngularVelocity;
using geometry::Displacement;
using geometry::Frame;
using geometry::InnerProduct;
using geometry::Permutation;
using geometry::Wedge;
using quantities::AngularFrequency;
using quantities::Speed;
using quantities::si::Day;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::VanishesBefore;

class RigidMotionTest : public testing::Test {
 protected:
  // Our Earth and Moon both run a bit slow.
  AngularFrequency const earth_rotation_speed_ = 2 * π * Radian / Day;
  AngularFrequency const moon_rotation_speed_ = 2 * π * Radian / (30 * Day);
  Length const earth_moon_distance_ = 384'400 * Kilo(Metre);
  // Nonrotating frame fixing the centre of the Earth.  The North pole is the
  // positive z axis, the x axis points towards the Moon,
  // the reference frame is right-handed.
  using Geocentric = Frame<serialization::Frame::TestTag,
                           serialization::Frame::TEST, true>;
  // Nonrotating frame fixing the centre of the Moon.  The North pole is the
  // positive z axis, the y axis points away from the Earth,
  // the reference frame is left-handed.
  using Selenocentric = Frame<serialization::Frame::TestTag,
                              serialization::Frame::TEST1, true>;
  // Rotating frame fixing the Earth's surface.  The North pole is the
  // positive z axis, the x axis points towards the Moon,
  // the reference frame is right-handed.
  using Terrestrial = Frame<serialization::Frame::TestTag,
                            serialization::Frame::TEST2, true>;
  // Rotating frame fixing the Moon's surface.
  // Nonrotating frame fixing the centre of the Moon.  The North pole is the
  // positive z axis, the y axis points away from the Earth,
  // the reference frame is left-handed.
  using Lunar = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST3, true>;

  AngularVelocity<Geocentric> const earth_rotation_ =
      AngularVelocity<Geocentric>(
          {0 * Radian / Second, 0 * Radian / Second, earth_rotation_speed_});
  RigidMotion<Geocentric, Terrestrial> const geocentric_to_terrestrial_ =
      RigidMotion<Geocentric, Terrestrial>(
          RigidTransformation<Geocentric, Terrestrial>(
              Geocentric::origin,
              Terrestrial::origin,
              OrthogonalMap<Geocentric, Terrestrial>::Identity()),
          earth_rotation_,
          Velocity<Geocentric>());

  // Our Moon is in an equatorial orbit for simplicity.
  AngularVelocity<Geocentric> const moon_orbit_ = AngularVelocity<Geocentric>(
      {0 * Radian / Second, 0 * Radian / Second, moon_rotation_speed_});
  Displacement<Geocentric> const earth_to_moon_ =
      Displacement<Geocentric>({earth_moon_distance_, 0 * Metre, 0 * Metre});

  RigidMotion<Geocentric, Selenocentric> const geocentric_to_selenocentric_ =
      RigidMotion<Geocentric, Selenocentric>(
          RigidTransformation<Geocentric, Selenocentric>(
              Geocentric::origin + earth_to_moon_,
              Selenocentric::origin,
              Permutation<Geocentric, Selenocentric>(
                  Permutation<Geocentric, Selenocentric>::YXZ).Forget()),
          AngularVelocity<Geocentric>(),
          moon_orbit_* earth_to_moon_ / Radian);

  AngularVelocity<Selenocentric> const moon_rotation_ =
      AngularVelocity<Selenocentric>(
          {0 * Radian / Second, 0 * Radian / Second, -moon_rotation_speed_});
  RigidMotion<Selenocentric, Lunar> const selenocentric_to_lunar_ =
      RigidMotion<Selenocentric, Lunar>(
          RigidTransformation<Selenocentric, Lunar>(
              Selenocentric::origin,
              Lunar::origin,
              OrthogonalMap<Selenocentric, Lunar>::Identity()),
          moon_rotation_,
          Velocity<Selenocentric>());

  // General degrees of freedom.
  DegreesOfFreedom<Terrestrial> const degrees_of_freedom_ = {
      Terrestrial::origin +
          Displacement<Terrestrial>({earth_moon_distance_ / 3,
                                     -earth_moon_distance_ / 5,
                                     3 * earth_moon_distance_ / 7}),
      earth_rotation_speed_ * earth_moon_distance_ / Radian *
          Vector<double, Terrestrial>({-0.5, 0.42, 2.1})};
};

TEST_F(RigidMotionTest, TidalLocking) {
  RigidMotion<Geocentric, Lunar> const geocentric_to_lunar =
      selenocentric_to_lunar_ * geocentric_to_selenocentric_;
  DegreesOfFreedom<Lunar> const earth_degrees_of_freedom =
      geocentric_to_lunar({Geocentric::origin, Velocity<Geocentric>()});
  EXPECT_EQ(Displacement<Lunar>({0 * Metre, -earth_moon_distance_, 0 * Metre}),
            earth_degrees_of_freedom.position() - Lunar::origin);
  Speed const moon_speed = (moon_orbit_ * earth_to_moon_ / Radian).Norm();
  EXPECT_THAT(earth_degrees_of_freedom.velocity(),
              Componentwise(VanishesBefore(moon_speed, 0),
                            VanishesBefore(moon_speed, 3),
                            0 * Metre / Second));
}

TEST_F(RigidMotionTest, ApparentMoon) {
  RigidMotion<Selenocentric, Terrestrial> const selenocentric_to_terrestrial =
      geocentric_to_terrestrial_ * geocentric_to_selenocentric_.Inverse();
  DegreesOfFreedom<Terrestrial> const moon_degrees_of_freedom =
      selenocentric_to_terrestrial(
          {Selenocentric::origin, Velocity<Selenocentric>()});
  Displacement<Terrestrial> const earth_to_moon =
      moon_degrees_of_freedom.position() - Terrestrial::origin;
  EXPECT_EQ(
      Displacement<Terrestrial>({earth_moon_distance_, 0 * Metre, 0 * Metre}),
      earth_to_moon);
  AngularVelocity<Terrestrial> const moon_angular_velocity =
      Wedge(earth_to_moon, moon_degrees_of_freedom.velocity()) /
      earth_to_moon.Norm²() * Radian;
  EXPECT_THAT(moon_angular_velocity / (2 * π * Radian) * Day,
              Componentwise(0, 0, AlmostEquals(-29.0 / 30.0, 5)));
}

TEST_F(RigidMotionTest, GroupoidAssociativity) {
  auto const terrestrial_to_geocentric = geocentric_to_terrestrial_.Inverse();
  DegreesOfFreedom<Lunar> const d1 =
      ((selenocentric_to_lunar_ * geocentric_to_selenocentric_) *
       terrestrial_to_geocentric)(degrees_of_freedom_);
  DegreesOfFreedom<Lunar> const d2 =
      (selenocentric_to_lunar_ *
       (geocentric_to_selenocentric_ * terrestrial_to_geocentric))(
          degrees_of_freedom_);
  EXPECT_THAT(d1.position() - Lunar::origin,
              AlmostEquals(d2.position() - Lunar::origin, 0));
  EXPECT_THAT(d1.velocity(), AlmostEquals(d2.velocity(), 4));
}

TEST_F(RigidMotionTest, GroupoidAction) {
  auto const terrestrial_to_geocentric = geocentric_to_terrestrial_.Inverse();
  auto const geocentric_to_lunar =
      selenocentric_to_lunar_ * geocentric_to_selenocentric_;
  DegreesOfFreedom<Lunar> const d1 =
      (geocentric_to_lunar * terrestrial_to_geocentric)(degrees_of_freedom_);
  DegreesOfFreedom<Lunar> const d2 =
      geocentric_to_lunar(terrestrial_to_geocentric(degrees_of_freedom_));
  EXPECT_THAT(d1.position() - Lunar::origin,
              AlmostEquals(d2.position() - Lunar::origin, 5));
  EXPECT_THAT(d1.velocity(), AlmostEquals(d2.velocity(), 1));
}

TEST_F(RigidMotionTest, GroupoidInverse) {
  auto const terrestrial_to_lunar = selenocentric_to_lunar_ *
                                    geocentric_to_selenocentric_ *
                                    geocentric_to_terrestrial_.Inverse();
  DegreesOfFreedom<Terrestrial> const d1 =
      (terrestrial_to_lunar.Inverse() *
       terrestrial_to_lunar)(degrees_of_freedom_);
  DegreesOfFreedom<Terrestrial> const d2 =
      terrestrial_to_lunar.Inverse()(terrestrial_to_lunar(degrees_of_freedom_));
  EXPECT_THAT(
      d1.position() - Terrestrial::origin,
      AlmostEquals(degrees_of_freedom_.position() - Terrestrial::origin, 0));
  EXPECT_THAT(d1.velocity(), AlmostEquals(degrees_of_freedom_.velocity(), 4));
  EXPECT_THAT(
      d2.position() - Terrestrial::origin,
      AlmostEquals(degrees_of_freedom_.position() - Terrestrial::origin, 4));
  EXPECT_THAT(d2.velocity(), AlmostEquals(degrees_of_freedom_.velocity(), 6));
}

TEST_F(RigidMotionTest, SecondConstructor) {
  DegreesOfFreedom<Selenocentric> const earth_dof_in_selenocentric =
      geocentric_to_selenocentric_({Geocentric::origin,
                                    Velocity<Geocentric>{}});
  AngularVelocity<Selenocentric> const earth_rotation_in_selenocentric;

  RigidMotion<Geocentric, Selenocentric> const geocentric_to_selenocentric(
      geocentric_to_selenocentric_.rigid_transformation(),
      earth_rotation_in_selenocentric,
      earth_dof_in_selenocentric.velocity());

  DegreesOfFreedom<Geocentric> const degrees_of_freedom = {
      Geocentric::origin +
          Displacement<Geocentric>({earth_moon_distance_ / 3,
                                    -earth_moon_distance_ / 5,
                                    3 * earth_moon_distance_ / 7}),
      earth_rotation_speed_ * earth_moon_distance_ / Radian *
          Vector<double, Geocentric>({-0.5, 0.42, 2.1})};

  auto const d1 = geocentric_to_selenocentric(degrees_of_freedom);
  auto const d2 = geocentric_to_selenocentric_(degrees_of_freedom);
  EXPECT_THAT(d1.position() - Selenocentric::origin,
              AlmostEquals(d2.position() - Selenocentric::origin, 0));
  EXPECT_THAT(d1.velocity(), AlmostEquals(d2.velocity(), 0));
}

}  // namespace internal_rigid_motion
}  // namespace physics
}  // namespace principia
