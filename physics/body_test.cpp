#include "physics/body.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using principia::geometry::Normalize;
using testing::IsNull;
using testing::NotNull;

namespace principia {
namespace physics {

class BodyTest : public testing::Test {
 protected:
  enum class Tag {
    kWorld,
  };

  using World = Frame<Tag, Tag::kWorld, true>;

  // We need that so the comma doesn't get caught in macros.
  using Direction = Vector<double, World>;

  Direction axis_ = Normalize(Direction({-1, 2, 5}));
  MasslessBody massless_body_;
  MassiveBody massive_body_ =
      MassiveBody(42 * SIUnit<GravitationalParameter>());
  OblateBody<World> oblate_body_ =
      OblateBody<World>(17 * SIUnit<GravitationalParameter>(),
                        163 * SIUnit<Order2ZonalCoefficient>(),
                        axis_);
};

using BodyDeathTest = BodyTest;

TEST_F(BodyTest, MasslessSerializationSuccess) {
  serialization::Body message;
  MasslessBody const* cast_massless_body;
  massless_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massless_body());
  EXPECT_FALSE(message.has_massive_body());

  // Direct deserialization.
  // No members to test in this class, we just check that it doesn't crash.
  massless_body_ = *MasslessBody::ReadFromMessage(message);

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  // NOTE(egg): The &* is a quick way to explicitly forget |not_null|ness. We
  // cannot strip the |not_null| from the previous line because MSVC does not
  // support move conversion at the moment.
  cast_massless_body = dynamic_cast<MasslessBody const*>(&*body);
  EXPECT_THAT(cast_massless_body, NotNull());
}

// The best serialization revenge.
TEST_F(BodyTest, MassiveSerializationSuccess) {
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
  cast_massive_body = dynamic_cast<MassiveBody*>(&*body);
  EXPECT_THAT(cast_massive_body, NotNull());
  EXPECT_EQ(massive_body_.gravitational_parameter(),
            cast_massive_body->gravitational_parameter());
}

TEST_F(BodyTest, OblateSerializationSuccess) {
  serialization::Body message;
  OblateBody<UnknownInertialFrame> const* unknown_cast_oblate_body;
  OblateBody<World> const* cast_oblate_body;
  oblate_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massive_body());
  EXPECT_FALSE(message.has_massless_body());
  EXPECT_TRUE(
      message.massive_body().HasExtension(
          serialization::OblateBody::oblate_body));
  EXPECT_EQ(17, message.massive_body().gravitational_parameter().magnitude());
  serialization::OblateBody const oblateness_information =
      message.massive_body().GetExtension(
          serialization::OblateBody::oblate_body);
  EXPECT_EQ(163, oblateness_information.j2().magnitude());
  EXPECT_EQ(axis_, Direction::ReadFromMessage(oblateness_information.axis()));

  // Direct deserialization.
  OblateBody<World> const oblate_body =
      *OblateBody<World>::ReadFromMessage(message);
  EXPECT_EQ(oblate_body_.gravitational_parameter(), 
            oblate_body.gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), oblate_body.j2());
  EXPECT_EQ(oblate_body_.axis(), oblate_body.axis());

  // Dispatching from |MassiveBody|.
  not_null<std::unique_ptr<MassiveBody const>> const massive_body =
      MassiveBody::ReadFromMessage(message);
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            massive_body->gravitational_parameter());
  unknown_cast_oblate_body =
      dynamic_cast<OblateBody<UnknownInertialFrame> const*>(&*massive_body);
  EXPECT_THAT(unknown_cast_oblate_body, NotNull());
  cast_oblate_body = dynamic_cast<OblateBody<World> const*>(&*massive_body);
  EXPECT_THAT(cast_oblate_body, IsNull());
  cast_oblate_body =
      reinterpret_cast<OblateBody<World> const*>(unknown_cast_oblate_body);
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            cast_oblate_body->gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), cast_oblate_body->j2());
  EXPECT_EQ(oblate_body_.axis(), cast_oblate_body->axis());

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  cast_oblate_body = dynamic_cast<OblateBody<World> const*>(&*body);
  unknown_cast_oblate_body =
      dynamic_cast<OblateBody<UnknownInertialFrame> const*>(&*body);
  EXPECT_THAT(unknown_cast_oblate_body, NotNull());
  cast_oblate_body = dynamic_cast<OblateBody<World> const*>(&*body);
  EXPECT_THAT(cast_oblate_body, IsNull());
  cast_oblate_body =
      reinterpret_cast<OblateBody<World> const*>(unknown_cast_oblate_body);
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            cast_oblate_body->gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), cast_oblate_body->j2());
  EXPECT_EQ(oblate_body_.axis(), cast_oblate_body->axis());;
}

}  // namespace physics
}  // namespace principia
