#include "base/pull_serializer.hpp"

#include <vector>

#include "gmock/gmock.h"
#include "serialization/physics.pb.h"

namespace principia {

using serialization::Pair;
using serialization::Point;
using serialization::Quantity;
using serialization::Trajectory;
using ::testing::ElementsAreArray;

namespace base {

class PullSerializerTest : public ::testing::Test {
 protected:
  int const kChunkSize = 99;

  PullSerializerTest()
      : pull_serializer_(kChunkSize) {
    // Build a biggish protobuf for serialization.
    for (int i = 0; i < 100; ++i) {
      Trajectory::InstantaneousDegreesOfFreedom* idof =
          trajectory_.add_timeline();
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
  }

  PullSerializer pull_serializer_;
  Trajectory trajectory_;
};

TEST_F(PullSerializerTest, Test) {
  pull_serializer_.Start(&trajectory_);
  std::vector<int> actual_sizes;
  std::vector<int> expected_sizes(53, kChunkSize);
  expected_sizes.push_back(53);
  for (;;) {
    PullSerializer::Data const data = pull_serializer_.Pull();
    if (data.size == 0) {
      break;
    }
    actual_sizes.push_back(data.size);
  }
  EXPECT_THAT(actual_sizes, ElementsAreArray(expected_sizes));
}

}  // namespace base
}  // namespace principia
