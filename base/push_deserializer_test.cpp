#include "base/push_deserializer.hpp"

#include <algorithm>
#include <cstring>
#include <list>
#include <string>
#include <thread>  // NOLINT(build/c++11)
#include <vector>

#include "base/bytes.hpp"
#include "base/not_null.hpp"
#include "base/pull_serializer.hpp"
#include "gmock/gmock.h"
#include "serialization/physics.pb.h"

namespace principia {

using serialization::Pair;
using serialization::Point;
using serialization::Quantity;
using serialization::Trajectory;
using ::testing::ElementsAreArray;

namespace base {

namespace {
int const kDeserializerChunkSize = 99;
int const kRunsPerTest = 1000;
int const kSerializerChunkSize = 99;
int const kNumberOfChunks = 3;
const char kStart[] = "START";
}  // namespace

class PushDeserializerTest : public ::testing::Test {
 protected:
  PushDeserializerTest()
      : pull_serializer_(std::make_unique<PullSerializer>(
                             kSerializerChunkSize, kNumberOfChunks)),
        push_deserializer_(std::make_unique<PushDeserializer>(
                               kDeserializerChunkSize, kNumberOfChunks)),
        stream_(std::bind(&PushDeserializerTest::OnEmpty, this, &strings_)) {}

  static not_null<std::unique_ptr<Trajectory const>> BuildTrajectory() {
    not_null<std::unique_ptr<Trajectory>> result =
        make_not_null_unique<Trajectory>();
    // Build a biggish protobuf for serialization.
    for (int i = 0; i < 100; ++i) {
      Trajectory::InstantaneousDegreesOfFreedom* idof = result->add_timeline();
      Point* instant = idof->mutable_instant();
      Quantity* scalar = instant->mutable_scalar();
      scalar->set_dimensions(3);
      scalar->set_magnitude(3 * i);
      Pair* dof = idof->mutable_degrees_of_freedom();
      Pair::Element* t1 = dof->mutable_t1();
      Point* point1 = t1->mutable_point();
      Quantity* scalar1 = point1->mutable_scalar();
      scalar1->set_dimensions(1);
      scalar1->set_magnitude(i);
      Pair::Element* t2 = dof->mutable_t2();
      Point* point2 = t2->mutable_point();
      Quantity* scalar2 = point2->mutable_scalar();
      scalar2->set_dimensions(2);
      scalar2->set_magnitude(2 * i);
    }
    return std::move(result);
  }

  static void Stomp(Bytes const& bytes) {
    std::memset(bytes.data, 0xCD, static_cast<size_t>(bytes.size));
  }

  // Returns the first string in the list.  Note that the very first string is
  // always discarded.
  Bytes OnEmpty(not_null<std::list<std::string>*> const strings) {
    strings->pop_front();
    CHECK(!strings->empty());
    std::string& front = strings->front();
    return Bytes(reinterpret_cast<std::uint8_t*>(&front[0]),
                 static_cast<std::int64_t>(front.size()));
  }

  std::unique_ptr<PullSerializer> pull_serializer_;
  std::unique_ptr<PushDeserializer> push_deserializer_;
  internal::DelegatingArrayInputStream stream_;
  std::list<std::string> strings_;
};

using PushDeserializerDeathTest = PushDeserializerTest;

TEST_F(PushDeserializerTest, Stream) {
  void const* data;
  int size;

  strings_ = {kStart, "abc"};
  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(3, size);
  EXPECT_EQ("abc", std::string(static_cast<char const*>(data), size));
  EXPECT_EQ(3, stream_.ByteCount());

  strings_ = {kStart, ""};
  EXPECT_FALSE(stream_.Next(&data, &size));
  EXPECT_EQ(3, stream_.ByteCount());

  strings_ = {kStart, "abc", "xyzt"};
  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(3, size);
  EXPECT_EQ("abc", std::string(static_cast<char const*>(data), size));
  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(4, size);
  EXPECT_EQ("xyzt", std::string(static_cast<char const*>(data), size));
  EXPECT_EQ(10, stream_.ByteCount());

  strings_ = {kStart, "abc", "xyzt", "uvw", ""};
  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(3, size);
  EXPECT_TRUE(stream_.Skip(2));
  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(2, size);
  EXPECT_EQ("zt", std::string(static_cast<char const*>(data), size));
  EXPECT_FALSE(stream_.Skip(5));
  EXPECT_EQ(20, stream_.ByteCount());

  strings_ = {kStart, "abc"};
  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(3, size);
  stream_.BackUp(1);
  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(1, size);
  EXPECT_EQ("c", std::string(static_cast<char const*>(data), size));
  EXPECT_EQ(23, stream_.ByteCount());
}

TEST_F(PushDeserializerTest, DeserializationThreading) {
  Trajectory read_trajectory;
  auto const trajectory = BuildTrajectory();
  int const byte_size = trajectory->ByteSize();
  auto serialized_trajectory =
      std::make_unique<std::uint8_t[]>(byte_size);

  for (int i = 0; i < kRunsPerTest; ++i) {
    push_deserializer_ = std::make_unique<PushDeserializer>(
        kDeserializerChunkSize, kNumberOfChunks);

    trajectory->SerializePartialToArray(&serialized_trajectory[0], byte_size);
    push_deserializer_->Start(&read_trajectory);
    Bytes bytes(serialized_trajectory.get(), byte_size);
    push_deserializer_->Push(bytes,
                             std::bind(&PushDeserializerTest::Stomp, bytes));
    push_deserializer_->Push(Bytes(), nullptr);

    // Destroying the deserializer waits until deserialization is done.
    push_deserializer_.reset();
  }
}

// Exercise concurrent serialization and deserialization.
TEST_F(PushDeserializerTest, SerializationDeserialization) {
  Trajectory read_trajectory;
  auto const trajectory = BuildTrajectory();
  int const byte_size = trajectory->ByteSize();
  for (int i = 0; i < kRunsPerTest; ++i) {
    auto trajectory = BuildTrajectory();
    auto storage = std::make_unique<std::uint8_t[]>(byte_size);
    std::uint8_t* data = &storage[0];

    pull_serializer_ =
        std::make_unique<PullSerializer>(kSerializerChunkSize, kNumberOfChunks);
    push_deserializer_ = std::make_unique<PushDeserializer>(
        kDeserializerChunkSize, kNumberOfChunks);

    pull_serializer_->Start(std::move(trajectory));
    push_deserializer_->Start(&read_trajectory);
    for (;;) {
      Bytes const bytes = pull_serializer_->Pull();
      std::memcpy(data, bytes.data, static_cast<size_t>(bytes.size));
      push_deserializer_->Push(Bytes(data, bytes.size),
                               std::bind(&PushDeserializerTest::Stomp,
                                         Bytes(data, bytes.size)));
      data = &data[bytes.size];
      if (bytes.size == 0) {
        break;
      }
    }

    // Destroying the deserializer waits until deserialization is done.  It is
    // important that this happens before |storage| is destroyed.
    pull_serializer_.reset();
    push_deserializer_.reset();
  }
}

// Check that deserialization fails if we stomp on one extra bytes.
TEST_F(PushDeserializerDeathTest, Stomp) {
  EXPECT_DEATH({
    const int kStompChunk = 77;
    Trajectory read_trajectory;
    auto const trajectory = BuildTrajectory();
    int const byte_size = trajectory->ByteSize();
    auto serialized_trajectory =
        std::make_unique<std::uint8_t[]>(byte_size);
    trajectory->SerializePartialToArray(&serialized_trajectory[0], byte_size);
    push_deserializer_->Start(&read_trajectory);
    int left = byte_size;
    for (int i = 0; i < byte_size; i += kStompChunk) {
      Bytes bytes(&serialized_trajectory[i], std::min(left, kStompChunk));
      push_deserializer_->Push(bytes,
                               std::bind(&PushDeserializerTest::Stomp,
                                         Bytes(bytes.data, bytes.size + 1)));
      left -= kStompChunk;
    }
    push_deserializer_->Push(Bytes(), nullptr);
    push_deserializer_.reset();
    }, "failed.*Parse");
}

}  // namespace base
}  // namespace principia
