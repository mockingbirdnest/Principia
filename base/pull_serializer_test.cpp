#include "base/pull_serializer.hpp"

#include <list>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "serialization/physics.pb.h"

namespace principia {

using serialization::Pair;
using serialization::Point;
using serialization::Quantity;
using serialization::Trajectory;
using ::std::placeholders::_1;
using ::testing::ElementsAre;
using ::testing::ElementsAreArray;

namespace base {

namespace {
  int const kChunkSize = 99;
  int const kNumberOfChunks = 3;
  int const kRunsPerTest = 1000;
  int const kSmallChunkSize = 3;
}  // namespace

class PullSerializerTest : public ::testing::Test {
 protected:
  PullSerializerTest()
      : pull_serializer_(
            std::make_unique<PullSerializer>(kChunkSize, kNumberOfChunks)),
        stream_(Bytes(data_, kSmallChunkSize),
                std::bind(&PullSerializerTest::OnFull, this, _1, &strings_)) {}

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

  // Returns the first string in the list.  Note that the very first string is
  // always discarded.
  Bytes OnFull(Bytes const bytes,
               not_null<std::list<std::string>*> const strings) {
    strings->push_back(
        std::string(reinterpret_cast<const char*>(&bytes.data[0]),
                    static_cast<size_t>(bytes.size)));
    return Bytes(data_, kSmallChunkSize);
  }

  std::unique_ptr<PullSerializer> pull_serializer_;
  internal::DelegatingArrayOutputStream stream_;
  std::list<std::string> strings_;
  std::uint8_t data_[kSmallChunkSize];
};

TEST_F(PullSerializerTest, Stream) {
  void* data;
  int size;

  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(3, size);
  EXPECT_EQ(3, stream_.ByteCount());
  memcpy(data, "abc", 3);
  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(3, size);
  EXPECT_EQ(6, stream_.ByteCount());
  memcpy(data, "xy", 2);
  stream_.BackUp(1);
  EXPECT_EQ(5, stream_.ByteCount());
  EXPECT_TRUE(stream_.Next(&data, &size));
  EXPECT_EQ(3, size);
  EXPECT_EQ(8, stream_.ByteCount());
  memcpy(data, "uvw", 3);
  stream_.BackUp(2);
  EXPECT_EQ(6, stream_.ByteCount());
  EXPECT_THAT(strings_, ElementsAre("abc", "xy", "u"));
}

TEST_F(PullSerializerTest, SerializationSizes) {
  auto trajectory = BuildTrajectory();
  pull_serializer_->Start(std::move(trajectory));
  std::vector<std::int64_t> actual_sizes;
  std::vector<std::int64_t> expected_sizes(53, kChunkSize);
  expected_sizes.push_back(53);
  for (;;) {
    Bytes const bytes = pull_serializer_->Pull();
    if (bytes.size == 0) {
      break;
    }
    actual_sizes.push_back(bytes.size);
  }
  EXPECT_THAT(actual_sizes, ElementsAreArray(expected_sizes));
}

TEST_F(PullSerializerTest, SerializationThreading) {
  Trajectory read_trajectory;
  auto const trajectory = BuildTrajectory();
  int const byte_size = trajectory->ByteSize();
  auto expected_serialized_trajectory =
      std::make_unique<std::uint8_t[]>(byte_size);
  trajectory->SerializePartialToArray(&expected_serialized_trajectory[0],
                                      byte_size);

  // Run this test repeatedly to detect threading issues (it will flake in case
  // of problems).
  for (int i = 0; i < kRunsPerTest; ++i) {
    auto trajectory = BuildTrajectory();
    auto actual_serialized_trajectory =
        std::make_unique<std::uint8_t[]>(byte_size);
    std::uint8_t* data = &actual_serialized_trajectory[0];

    // The serialization happens concurrently with the test.
    pull_serializer_ =
        std::make_unique<PullSerializer>(kChunkSize, kNumberOfChunks);
    pull_serializer_->Start(std::move(trajectory));
    for (;;) {
      Bytes const bytes = pull_serializer_->Pull();
      std::memcpy(data, bytes.data, static_cast<size_t>(bytes.size));
      data = &data[bytes.size];
      if (bytes.size == 0) {
        break;
      }
    }
    pull_serializer_.reset();

    // Check if the serialized version can be parsed and if not print the first
    // difference.
    if (!read_trajectory.ParseFromArray(&actual_serialized_trajectory[0],
                                        byte_size)) {
      for (int i = 0; i < byte_size; ++i) {
        if (expected_serialized_trajectory[i] !=
            actual_serialized_trajectory[i]) {
          LOG(FATAL) << "position=" << i
                     << ", expected="
                     << static_cast<int>(expected_serialized_trajectory[i])
                     << ", actual="
                     << static_cast<int>(actual_serialized_trajectory[i]);
        }
      }
    }
  }
}

}  // namespace base
}  // namespace principia
