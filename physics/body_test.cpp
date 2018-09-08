
#include "physics/body.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace physics {
namespace internal_body {

using geometry::AngularVelocity;
using geometry::Frame;
using geometry::Instant;
using geometry::Normalize;
using geometry::RadiusLatitudeLongitude;
using geometry::Vector;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Degree2SphericalHarmonicCoefficient;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::SIUnit;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using ::testing::IsNull;
using ::testing::NotNull;

class BodyTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;

  // We need that so the comma doesn't get caught in macros.
  using Direction = Vector<double, World>;

  template<typename Tag, Tag tag>
  void TestRotatingBody() {
    using F = Frame<Tag, tag, true>;

    auto const rotating_body =
        RotatingBody<F>(17 * SIUnit<GravitationalParameter>(),
                        typename RotatingBody<F>::Parameters(
                            2 * Metre,
                            3 * Radian,
                            Instant() + 4 * Second,
                            angular_frequency_,
                            right_ascension_of_pole_,
                            declination_of_pole_));

    serialization::Body message;
    RotatingBody<F> const* cast_rotating_body;
    rotating_body.WriteToMessage(&message);
    EXPECT_TRUE(message.has_massive_body());
    EXPECT_FALSE(message.has_massless_body());
    EXPECT_TRUE(message.massive_body().HasExtension(
                    serialization::RotatingBody::extension));

    not_null<std::unique_ptr<MassiveBody const>> const massive_body =
        MassiveBody::ReadFromMessage(message);
    EXPECT_EQ(rotating_body.gravitational_parameter(),
              massive_body->gravitational_parameter());
    cast_rotating_body =
        dynamic_cast_not_null<RotatingBody<F> const*>(massive_body.get());
    EXPECT_THAT(cast_rotating_body, NotNull());
  }

  AngularFrequency const angular_frequency_ = -1.5 * Radian / Second;
  Angle const right_ascension_of_pole_ = 37 * Degree;
  Angle const declination_of_pole_ = 123 * Degree;

  MasslessBody massless_body_;
  MassiveBody massive_body_ =
      MassiveBody(42 * SIUnit<GravitationalParameter>());
  RotatingBody<World> rotating_body_ =
      RotatingBody<World>(17 * SIUnit<GravitationalParameter>(),
                          RotatingBody<World>::Parameters(
                              1 * Metre,
                              3 * Radian,
                              Instant() + 4 * Second,
                              angular_frequency_,
                              right_ascension_of_pole_,
                              declination_of_pole_));
  OblateBody<World> oblate_body_ =
      OblateBody<World>(17 * SIUnit<GravitationalParameter>(),
                        RotatingBody<World>::Parameters(
                            1 * Metre,
                            3 * Radian,
                            Instant() + 4 * Second,
                            angular_frequency_,
                            right_ascension_of_pole_,
                            declination_of_pole_),
                        OblateBody<World>::Parameters(
                            6,
                            1 * Metre));
};

using BodyDeathTest = BodyTest;

TEST_F(BodyTest, MasslessSerializationSuccess) {
  EXPECT_TRUE(massless_body_.is_massless());
  EXPECT_FALSE(massless_body_.is_oblate());

  serialization::Body message;
  MasslessBody const* cast_massless_body;
  massless_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massless_body());
  EXPECT_FALSE(message.has_massive_body());

  // Direct deserialization.
  // No members to test in this class, we just check that it doesn't crash.
  massless_body_ = *MasslessBody::ReadFromMessage(message);

  // Dispatching from |Body|.  Need two steps to add const and remove
  // |not_null|.
  not_null<std::unique_ptr<Body const>> body = Body::ReadFromMessage(message);
  cast_massless_body = dynamic_cast_not_null<MasslessBody const*>(body.get());
  EXPECT_THAT(cast_massless_body, NotNull());
}

// The best serialization revenge.
TEST_F(BodyTest, MassiveSerializationSuccess) {
  EXPECT_FALSE(massive_body_.is_massless());
  EXPECT_FALSE(massive_body_.is_oblate());

  serialization::Body message;
  MassiveBody const* cast_massive_body;
  massive_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massive_body());
  EXPECT_FALSE(message.has_massless_body());
  EXPECT_EQ(42, message.massive_body().gravitational_parameter().magnitude());

  // Direct deserialization.
  MassiveBody const massive_body = *MassiveBody::ReadFromMessage(message);
  EXPECT_EQ(massive_body_.gravitational_parameter(),
            massive_body.gravitational_parameter());

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body>> body = Body::ReadFromMessage(message);
  cast_massive_body = dynamic_cast_not_null<MassiveBody*>(body.get());
  EXPECT_THAT(cast_massive_body, NotNull());
  EXPECT_EQ(massive_body_.gravitational_parameter(),
            cast_massive_body->gravitational_parameter());
}

TEST_F(BodyTest, RotatingSerializationSuccess) {
  EXPECT_FALSE(rotating_body_.is_massless());
  EXPECT_FALSE(rotating_body_.is_oblate());

  serialization::Body message;
  RotatingBody<World> const* cast_rotating_body;
  rotating_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massive_body());
  EXPECT_FALSE(message.has_massless_body());
  EXPECT_TRUE(message.massive_body().HasExtension(
                  serialization::RotatingBody::extension));
  EXPECT_EQ(17, message.massive_body().gravitational_parameter().magnitude());
  serialization::RotatingBody const rotating_body_extension =
      message.massive_body().GetExtension(
          serialization::RotatingBody::extension);
  EXPECT_EQ(3, rotating_body_extension.reference_angle().magnitude());
  EXPECT_EQ(4,
            rotating_body_extension.reference_instant().scalar().magnitude());
  EXPECT_EQ(angular_frequency_,
            AngularFrequency::ReadFromMessage(
                rotating_body_extension.angular_frequency()));
  EXPECT_EQ(right_ascension_of_pole_,
            Angle::ReadFromMessage(
                rotating_body_extension.right_ascension_of_pole()));
  EXPECT_EQ(declination_of_pole_,
            Angle::ReadFromMessage(
                rotating_body_extension.declination_of_pole()));

  // Dispatching from |MassiveBody|.
  not_null<std::unique_ptr<MassiveBody const>> const massive_body =
      MassiveBody::ReadFromMessage(message);
  EXPECT_EQ(rotating_body_.gravitational_parameter(),
            massive_body->gravitational_parameter());
  cast_rotating_body =
      dynamic_cast_not_null<RotatingBody<World> const*>(massive_body.get());
  EXPECT_THAT(cast_rotating_body, NotNull());
  EXPECT_EQ(rotating_body_.gravitational_parameter(),
            cast_rotating_body->gravitational_parameter());
  EXPECT_EQ(rotating_body_.angular_velocity(),
            cast_rotating_body->angular_velocity());
  EXPECT_EQ(rotating_body_.AngleAt(Instant()),
            cast_rotating_body->AngleAt(Instant()));

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  cast_rotating_body =
      dynamic_cast_not_null<RotatingBody<World> const*>(body.get());
  EXPECT_THAT(cast_rotating_body, NotNull());
  EXPECT_EQ(rotating_body_.gravitational_parameter(),
            cast_rotating_body->gravitational_parameter());
  EXPECT_EQ(rotating_body_.angular_velocity(),
            cast_rotating_body->angular_velocity());
  EXPECT_EQ(rotating_body_.AngleAt(Instant()),
            cast_rotating_body->AngleAt(Instant()));
}

TEST_F(BodyTest, OblateSerializationSuccess) {
  EXPECT_FALSE(oblate_body_.is_massless());
  EXPECT_TRUE(oblate_body_.is_oblate());

  serialization::Body message;
  OblateBody<World> const* cast_oblate_body;
  oblate_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massive_body());
  EXPECT_FALSE(message.has_massless_body());
  EXPECT_TRUE(message.massive_body().GetExtension(
                  serialization::RotatingBody::extension).
                      HasExtension(serialization::OblateBody::extension));
  EXPECT_EQ(17, message.massive_body().gravitational_parameter().magnitude());
  serialization::OblateBody const oblate_body_extension =
      message.massive_body().GetExtension(
                  serialization::RotatingBody::extension).
                      GetExtension(serialization::OblateBody::extension);
  EXPECT_EQ(6, oblate_body_extension.j2());

  // Dispatching from |MassiveBody|.
  not_null<std::unique_ptr<MassiveBody const>> const massive_body =
      MassiveBody::ReadFromMessage(message);
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            massive_body->gravitational_parameter());
  cast_oblate_body =
      dynamic_cast_not_null<OblateBody<World> const*>(massive_body.get());
  EXPECT_THAT(cast_oblate_body, NotNull());
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            cast_oblate_body->gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), cast_oblate_body->j2());
  EXPECT_EQ(oblate_body_.polar_axis(), cast_oblate_body->polar_axis());

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  cast_oblate_body =
      dynamic_cast_not_null<OblateBody<World> const*>(body.get());
  EXPECT_THAT(cast_oblate_body, NotNull());
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            cast_oblate_body->gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), cast_oblate_body->j2());
  EXPECT_EQ(oblate_body_.polar_axis(), cast_oblate_body->polar_axis());
}

TEST_F(BodyTest, OblateSerializationCompatibility) {
  EXPECT_FALSE(oblate_body_.is_massless());
  EXPECT_TRUE(oblate_body_.is_oblate());

  // Construct a pre-Descartes message.
  serialization::Body message;
  OblateBody<World> const* cast_oblate_body;
  oblate_body_.WriteToMessage(&message);
  serialization::OblateBody* const oblate_body_extension =
      message.mutable_massive_body()->MutableExtension(
                  serialization::RotatingBody::extension)->
                      MutableExtension(serialization::OblateBody::extension);
  oblate_body_extension->clear_reference_radius();
  Degree2SphericalHarmonicCoefficient pre_descartes_j2 =
      7 * SIUnit<Degree2SphericalHarmonicCoefficient>();
  pre_descartes_j2.WriteToMessage(
      oblate_body_extension->mutable_pre_descartes_j2());

  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  cast_oblate_body =
      dynamic_cast_not_null<OblateBody<World> const*>(body.get());
  Length const reference_radius = 1 * Metre;
  EXPECT_EQ(7 / (cast_oblate_body->gravitational_parameter() /
                 SIUnit<GravitationalParameter>()),
            cast_oblate_body->j2());
  EXPECT_EQ(reference_radius, cast_oblate_body->reference_radius());
  EXPECT_EQ(7 * SIUnit<Degree2SphericalHarmonicCoefficient>() /
                cast_oblate_body->gravitational_parameter(),
            cast_oblate_body->j2_over_μ());
}

TEST_F(BodyTest, AllFrames) {
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::ALICE_SUN>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::ALICE_WORLD>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::BARYCENTRIC>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::NAVIGATION>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::WORLD>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::WORLD_SUN>();

  TestRotatingBody<serialization::Frame::SolarSystemTag,
                   serialization::Frame::GCRS>();
  TestRotatingBody<serialization::Frame::SolarSystemTag,
                   serialization::Frame::ICRS>();
  TestRotatingBody<serialization::Frame::SolarSystemTag,
                   serialization::Frame::ITRS>();

  TestRotatingBody<serialization::Frame::TestTag, serialization::Frame::TEST>();
  TestRotatingBody<serialization::Frame::TestTag,
                   serialization::Frame::TEST1>();
  TestRotatingBody<serialization::Frame::TestTag,
                   serialization::Frame::TEST2>();
  TestRotatingBody<serialization::Frame::TestTag, serialization::Frame::FROM>();
  TestRotatingBody<serialization::Frame::TestTag,
                   serialization::Frame::THROUGH>();
  TestRotatingBody<serialization::Frame::TestTag, serialization::Frame::TO>();
}

}  // namespace internal_body
}  // namespace physics
}  // namespace principia
