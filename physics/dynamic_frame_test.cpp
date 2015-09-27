#include "physics/dynamic_frame.hpp"

#include "geometry/frame.hpp"
#include "geometry/permutation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::AngularVelocity;
using geometry::Frame;
using geometry::Permutation;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::si::Day;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

namespace physics {

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
  RigidMotion<Terrestrial, Geocentric> const terrestrial_to_geocentric_ =
      RigidMotion<Terrestrial, Geocentric>(
          RigidTransformation<Terrestrial, Geocentric>(
              Terrestrial::origin,
              Geocentric::origin,
              OrthogonalMap<Terrestrial, Geocentric>::Identity()),
          earth_rotation_,
          Velocity<Geocentric>());

   // Our Moon is in an equatorial orbit for simplicity.
   AngularVelocity<Geocentric> const moon_orbit_ =
      AngularVelocity<Geocentric>(
          {0 * Radian / Second, 0 * Radian / Second, moon_rotation_speed_});
   Displacement<Geocentric> const earth_to_moon_ =
       Displacement<Geocentric>({earth_moon_distance_, 0 * Metre, 0 * Metre});

   RigidMotion<Selenocentric, Geocentric> const selenocentric_to_geocentric_ =
       RigidMotion<Selenocentric, Geocentric>(
           RigidTransformation<Selenocentric, Geocentric>(
               Selenocentric::origin,
               Geocentric::origin + earth_to_moon_,
               Permutation<Selenocentric, Geocentric>(
                   Permutation<Selenocentric, Geocentric>::YXZ)
                   .Forget()),
           AngularVelocity<Geocentric>(),
           moon_orbit_ * earth_to_moon_ / Radian);

  AngularVelocity<Selenocentric> const moon_rotation_ =
      AngularVelocity<Selenocentric>(
          {0 * Radian / Second, 0 * Radian / Second, -moon_rotation_speed_});
  RigidMotion<Lunar, Selenocentric> const lunar_to_selenocentric_ =
      RigidMotion<Lunar, Selenocentric>(
          RigidTransformation<Lunar, Selenocentric>(
              Lunar::origin,
              Selenocentric::origin,
              OrthogonalMap<Lunar, Selenocentric>::Identity()),
          moon_rotation_,
          Velocity<Selenocentric>());
};

TEST_F(RigidMotionTest, TidalLocking) {
  RigidMotion<Geocentric, Lunar> const geocentric_to_lunar =
      lunar_to_selenocentric_.Inverse() *
      selenocentric_to_geocentric_.Inverse();
  DegreesOfFreedom<Lunar> const earth_degrees_of_freedom =
      geocentric_to_lunar({Geocentric::origin, Velocity<Geocentric>()});
  EXPECT_EQ(earth_degrees_of_freedom.position() - Lunar::origin,
            Displacement<Lunar>({0 * Metre, -earth_moon_distance_, 0 * Metre}));
  EXPECT_EQ(earth_degrees_of_freedom.velocity(), Velocity<Lunar>());
}

}  // namespace physics
}  // namespace principia
