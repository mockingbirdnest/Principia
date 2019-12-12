
#include "geometry/frame.hpp"

#include "base/traits.hpp"
#include "glog/logging.h"
#include "google/protobuf/descriptor.h"
#include "gtest/gtest.h"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

class FrameTest : public testing::Test {
 protected:
  using World1 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST1,
                       Inertial,
                       Handedness::Right>;
  using World2 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST2,
                       Inertial,
                       Handedness::Right>;
  using World3 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST1,
                       NonInertial,
                       Handedness::Right>;
  using World4 = Frame<serialization::Frame::SolarSystemTag,
                       serialization::Frame::ICRS,
                       Inertial,
                       Handedness::Right>;

  using F1 = Frame<enum class F1Tag>;
  using F2 = Frame<enum class F2Tag>;
  using F3 = Frame<enum class F3Tag, F3Tag{}, Inertial>;
  static_assert(!std::is_same_v<F1, F2>);
  static_assert(!std::is_same_v<F1, F3>);
  static_assert(!std::is_same_v<F2, F3>);

  static_assert(base::is_serializable_v<World1>);
  static_assert(!base::is_serializable_v<F1>);
};

using FrameDeathTest = FrameTest;

// Uncomment to check that non-serializable frames are detected at compile-time.
#if 0
TEST_F(FrameTest, SerializationCompilationError) {
  serialization::Frame message;
  F1::WriteToMessage(&message);
  F2::WriteToMessage(&message);
  F3::WriteToMessage(&message);
}
#endif

TEST_F(FrameDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Frame message;
    World1::WriteToMessage(&message);
    World2::ReadFromMessage(message);
  }, R"(\(tag\) ==)");
  EXPECT_DEATH({
    serialization::Frame message;
    World1::WriteToMessage(&message);
    World3::ReadFromMessage(message);
  }, "is_inertial ==");
  EXPECT_DEATH({
    serialization::Frame message;
    World1::WriteToMessage(&message);
    World4::ReadFromMessage(message);
  }, "Fingerprint");
  EXPECT_DEATH({
    serialization::Frame message;
    World1::WriteToMessage(&message);
    message.set_tag_type_fingerprint(0xDEADBEEF);
    google::protobuf::EnumValueDescriptor const* enum_value_descriptor;
    bool is_inertial;
    ReadFrameFromMessage(message, enum_value_descriptor, is_inertial);
  }, "enum_value_descriptor");
  EXPECT_DEATH({
    serialization::Frame message;
    World1::WriteToMessage(&message);
    message.set_tag(666);
    google::protobuf::EnumValueDescriptor const* enum_value_descriptor;
    bool is_inertial;
    ReadFrameFromMessage(message, enum_value_descriptor, is_inertial);
  }, "enum_value_descriptor");
}

TEST_F(FrameTest, SerializationSuccess) {
  serialization::Frame message;
  World1::WriteToMessage(&message);
  EXPECT_TRUE(message.has_tag_type_fingerprint());
  // Be very worried if the following test fails: you have changed the
  // fingerprints and all your frames are belong to us.
  EXPECT_EQ(3211394049, message.tag_type_fingerprint());
  EXPECT_TRUE(message.has_tag());
  EXPECT_EQ(2, message.tag());
  EXPECT_TRUE(message.is_inertial());

  World1::ReadFromMessage(message);

  google::protobuf::EnumValueDescriptor const* enum_value_descriptor;
  bool is_inertial;
  ReadFrameFromMessage(message, enum_value_descriptor, is_inertial);
  EXPECT_EQ("principia.serialization.Frame.TEST1",
            enum_value_descriptor->full_name());
  EXPECT_TRUE(is_inertial);
}

}  // namespace geometry
}  // namespace principia
