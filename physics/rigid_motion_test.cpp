
#include "physics/rigid_motion.hpp"

#include "geometry/frame.hpp"
#include "geometry/permutation.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "geometry/signature.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {

using geometry::AngularVelocity;
using geometry::Arbitrary;
using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::OddPermutation;
using geometry::OrthogonalMap;
using geometry::Permutation;
using geometry::Point;
using geometry::Position;
using geometry::Quaternion;
using geometry::Rotation;
using geometry::Sign;
using geometry::Signature;
using geometry::Vector;
using geometry::Velocity;
using geometry::Wedge;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Speed;
using quantities::si::Day;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::IsNear;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::VanishesBefore;
using testing_utilities::operator""_⑴;

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
                           Inertial,
                           Handedness::Right,
                           serialization::Frame::TEST>;
  // Nonrotating frame fixing the centre of the Moon.  The North pole is the
  // positive z axis, the y axis points away from the Earth,
  // the reference frame is left-handed.
  using Selenocentric = Frame<serialization::Frame::TestTag,
                              Inertial,
                              Handedness::Left,
                              serialization::Frame::TEST1>;
  // Rotating frame fixing the Earth's surface.  The North pole is the
  // positive z axis, the x axis points towards the Moon,
  // the reference frame is right-handed.
  using Terrestrial = Frame<serialization::Frame::TestTag,
                            Arbitrary,
                            Handedness::Right,
                            serialization::Frame::TEST2>;
  // Rotating frame fixing the Moon's surface.
  // Nonrotating frame fixing the centre of the Moon.  The North pole is the
  // positive z axis, the y axis points away from the Earth,
  // the reference frame is left-handed.
  using Lunar = Frame<serialization::Frame::TestTag,
                      Arbitrary,
                      Handedness::Left,
                      serialization::Frame::TEST3>;

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
          Geocentric::unmoving);

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
              Permutation<Geocentric, Selenocentric>(OddPermutation::YXZ)
                  .Forget<OrthogonalMap>()),
          Geocentric::nonrotating,
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
          Selenocentric::unmoving);

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
      geocentric_to_lunar({Geocentric::origin, Geocentric::unmoving});
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
          {Selenocentric::origin, Selenocentric::unmoving});
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
  RigidMotion<Terrestrial, Selenocentric> const terrestrial_to_selenocentric1 =
      geocentric_to_selenocentric_ * geocentric_to_terrestrial_.Inverse();
  DegreesOfFreedom<Selenocentric> const terrestrial_dof_in_selenocentric =
      terrestrial_to_selenocentric1({Terrestrial::origin,
                                     Terrestrial::unmoving});
  AngularVelocity<Selenocentric> const terrestrial_rotation_in_selenocentric =
      terrestrial_to_selenocentric1.angular_velocity_of<Terrestrial>();

  RigidMotion<Terrestrial, Selenocentric> const terrestrial_to_selenocentric2(
      terrestrial_to_selenocentric1.rigid_transformation(),
      terrestrial_rotation_in_selenocentric,
      terrestrial_dof_in_selenocentric.velocity());

  auto const d1 = terrestrial_to_selenocentric1(degrees_of_freedom_);
  auto const d2 = terrestrial_to_selenocentric2(degrees_of_freedom_);
  EXPECT_THAT(d1.position() - Selenocentric::origin,
              AlmostEquals(d2.position() - Selenocentric::origin, 0));
  EXPECT_THAT(d1.velocity(), AlmostEquals(d2.velocity(), 3));
}

TEST_F(RigidMotionTest, Serialization) {
  serialization::RigidMotion message;

  auto const terrestrial_to_selenocentric1 =
      geocentric_to_selenocentric_ * geocentric_to_terrestrial_.Inverse();
  terrestrial_to_selenocentric1.WriteToMessage(&message);

  auto const terrestrial_to_selenocentric2 =
      RigidMotion<Terrestrial, Selenocentric>::ReadFromMessage(message);

  auto const d1 = terrestrial_to_selenocentric1(degrees_of_freedom_);
  auto const d2 = terrestrial_to_selenocentric2(degrees_of_freedom_);
  EXPECT_THAT(d1.position() - Selenocentric::origin,
              AlmostEquals(d2.position() - Selenocentric::origin, 0));
  EXPECT_THAT(d1.velocity(), AlmostEquals(d2.velocity(), 0));
}

TEST_F(RigidMotionTest, QuaternionNormalization) {
  using Barycentric =
      Frame<enum class BarycentricTag, Inertial, Handedness::Right>;
  using RigidPart =
      Frame<enum class RigidPartTag, Arbitrary, Handedness::Left>;
  using AliceWorld =
      Frame<enum class AliceWorldTag, Inertial, Handedness::Right>;
  using World = Frame<enum class WorldTag, Inertial, Handedness::Left>;

  AngularVelocity<RigidPart> const ω1(
      {-4.31524874936563274e-04 * Radian / Second,
       +3.58760764995860824e-03 * Radian / Second,
       -5.60639012817497860e-04 * Radian / Second});
  Velocity<RigidPart> const v1({+9.45391439618678553e-06 * Metre / Second,
                                -3.92345736956409459e-03 * Metre / Second,
                                -2.74220213842711266e-06 * Metre / Second});
  Position<RigidPart> const from1 = RigidPart::origin;
  Position<World> const to1 =
      World::origin + Displacement<World>({+3.37126515805721283e-02 * Metre,
                                           +1.21778284665197134e-03 * Metre,
                                           -2.06734146922826767e-02 * Metre});

  // Unity gives us quaternions that are normalized in single precision.  If we
  // don't normalize them in double precision, the tests below fail with huge
  // errors (orders of magnitude).
  Rotation<RigidPart, World> const rotation1(
      Normalize(Quaternion(-3.23732018470764160e-01,
                           {3.28581303358078003e-01,
                            -6.28009378910064697e-01,
                            -6.26766622066497803e-01})));

  RigidMotion<RigidPart, World> const rigid_motion1(
      RigidTransformation<RigidPart, World>(
          from1, to1, rotation1.Forget<OrthogonalMap>()),
      ω1,
      v1);

  AngularVelocity<World> const ω2({-1.09321577613522688e-36 * Radian / Second,
                                   +2.91570900559802080e-04 * Radian / Second,
                                   +5.42101086242752217e-20 * Radian / Second});
  Velocity<World> const v2({+3.10203149761506589e+06 * Metre / Second,
                            +2.41841229267265589e-02 * Metre / Second,
                            +2.45903920001263265e+06 * Metre / Second});
  Position<World> const from2 =
      World::origin + Displacement<World>({+4.91525322355270386e+05 * Metre,
                                           +1.01811877291461769e+03 * Metre,
                                           -3.44263258773803711e+05 * Metre});
  Position<Barycentric> const to2 =
      Barycentric::origin +
      Displacement<Barycentric>({-1.36081879406308479e+10 * Metre,
                                 +7.17490648244777508e+06 * Metre,
                                 -4.33826960235058723e+04 * Metre});
  Rotation<AliceWorld, Barycentric> const rotation2(
      Quaternion(6.36714825392272976e-01,
                 {-6.36714825392272976e-01,
                  -3.07561751727498445e-01,
                  +3.07561751727498445e-01}));
  OrthogonalMap<World, Barycentric> const orthogonal2 =
      rotation2.Forget<OrthogonalMap>() *
      Signature<World, AliceWorld>::CentralInversion().Forget<OrthogonalMap>();

  RigidMotion<World, Barycentric> const rigid_motion2(
      RigidTransformation<World, Barycentric>(from2, to2, orthogonal2), ω2, v2);

  DegreesOfFreedom<World> const d1 =
      rigid_motion1({RigidPart::origin, RigidPart::unmoving});
  DegreesOfFreedom<Barycentric> const d2 =
      (rigid_motion2 * rigid_motion1)({RigidPart::origin, RigidPart::unmoving});
  DegreesOfFreedom<Barycentric> const d3 =
      rigid_motion2(rigid_motion1({RigidPart::origin, RigidPart::unmoving}));
  DegreesOfFreedom<World> const d4 =
      rigid_motion2.Inverse()((rigid_motion2 * rigid_motion1)(
          {RigidPart::origin, RigidPart::unmoving}));
  DegreesOfFreedom<World> const d5 = rigid_motion2.Inverse()(
      rigid_motion2(rigid_motion1({RigidPart::origin, RigidPart::unmoving})));

  EXPECT_THAT(d2.position() - Barycentric::origin,
              AlmostEquals(d3.position() - Barycentric::origin, 0));
  EXPECT_THAT(d2.velocity(),
              RelativeErrorFrom(d3.velocity(), IsNear(1.1e-12_⑴)));
  EXPECT_THAT(d4.position() - World::origin,
              RelativeErrorFrom(d1.position() - World::origin,
                                IsNear(5.8e-6_⑴)));
  EXPECT_THAT(d4.velocity(),
              RelativeErrorFrom(d1.velocity(), IsNear(2.7e-6_⑴)));
  EXPECT_THAT(d5.position() - World::origin,
              RelativeErrorFrom(d1.position() - World::origin,
                                IsNear(5.8e-6_⑴)));
  EXPECT_THAT(d5.velocity(),
              RelativeErrorFrom(d1.velocity(), IsNear(5.7e-8_⑴)));
}

}  // namespace physics
}  // namespace principia
