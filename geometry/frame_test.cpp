#include "geometry/frame.hpp"

#include "base/concepts.hpp"
#include "glog/logging.h"
#include "google/protobuf/descriptor.h"
#include "gtest/gtest.h"
#include "serialization/geometry.pb.h"
#include "testing_utilities/check_well_formedness.hpp"  // 🧙 For PRINCIPIA_CHECK_ILL_FORMED.

namespace principia {
namespace geometry {

using namespace principia::base::_concepts;
using namespace principia::geometry::_frame;

class FrameTest : public testing::Test {
 public:
  using World1 = Frame<serialization::Frame::TestTag,
                       Inertial,
                       Handedness::Right,
                       serialization::Frame::TEST1>;
  using World2 = Frame<serialization::Frame::TestTag,
                       Inertial,
                       Handedness::Right,
                       serialization::Frame::TEST2>;
  using World3 = Frame<serialization::Frame::TestTag,
                       Arbitrary,
                       Handedness::Right,
                       serialization::Frame::TEST1>;
  using World4 = Frame<serialization::Frame::SolarSystemTag,
                       Inertial,
                       Handedness::Right,
                       serialization::Frame::ICRS>;

  using F1 = Frame<struct F1Tag>;
  using F2 = Frame<struct F2Tag>;
  using F3 = Frame<struct F3Tag, Inertial>;
  static_assert(!std::is_same_v<F1, F2>);
  static_assert(!std::is_same_v<F1, F3>);
  static_assert(!std::is_same_v<F2, F3>);

  static_assert(serializable<World1>);
  static_assert(!serializable<F1>);
};

using FrameDeathTest = FrameTest;

// Check that non-serializable frames are detected at compile-time.
PRINCIPIA_CHECK_WELL_FORMED_WITH_TYPES(
    World1::ReadFromMessage(message),
    (typename World1 = FrameTest::World1),
    with_variable<serialization::Frame> message);
PRINCIPIA_CHECK_ILL_FORMED_WITH_TYPES(
    F1::ReadFromMessage(message),
    (typename F1 = FrameTest::F1),
    with_variable<serialization::Frame> message);
PRINCIPIA_CHECK_ILL_FORMED_WITH_TYPES(
    F2::ReadFromMessage(message),
    (typename F2 = FrameTest::F2),
    with_variable<serialization::Frame> message);
PRINCIPIA_CHECK_ILL_FORMED_WITH_TYPES(
    F3::ReadFromMessage(message),
    (typename F3 = FrameTest::F3),
    with_variable<serialization::Frame> message);

TEST_F(FrameDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Frame message;
    World1::WriteToMessage(&message);
    World2::ReadFromMessage(message);
  }, R"(\(tag\(\)\) ==)");
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
