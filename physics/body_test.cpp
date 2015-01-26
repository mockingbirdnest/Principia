#include "physics/body.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using principia::geometry::Normalize;

namespace principia {
namespace physics {

class BodyTest : public testing::Test {
 protected:
  enum class Tag {
    kWorld,
  };

  using World = Frame<Tag, Tag::kWorld, true>;

  Vector<double, World> axis_ = Normalize(Vector<double, World>({-1, 2, 5}));
  MasslessBody massless_body_;
  MassiveBody massive_body_ =
      MassiveBody(42 * SIUnit<GravitationalParameter>());
  OblateBody<World> oblate_body_ =
      OblateBody<World>(17 * SIUnit<GravitationalParameter>(),
                        163 * SIUnit<Order2ZonalCoefficient>(),
                        axis_);
};

TEST_F(BodyTest, MasslessSerializationSuccess) {
  serialization::Body message;
  massless_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massless_body());
  EXPECT_FALSE(message.has_massive_body());
  massless_body_ = *MasslessBody::ReadFromMessage(message);
  not_null<std::unique_ptr<Body>> body = Body::ReadFromMessage(message);
  // NOTE(egg): The &* forgets |not_null|ness. We can't strip the |not_null|
  // from the previous line because MSVC does not support move conversion at the
  // moment.
  massless_body_ = *dynamic_cast<MasslessBody*>(&*body);
}

using BodyDeathTest = BodyTest;

}  // namespace physics
}  // namespace principia
