#include "base/serializer.hpp"

#include "gmock/gmock.h"
#include "serialization/physics.pb.h"

namespace principia {

using serialization::Pair;
using serialization::Point;
using serialization::Quantity;
using serialization::Trajectory;

namespace base {

class SerializerTest : public ::testing::Test {
 protected:
  SerializerTest()
      : serializer_(12) {
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

  Serializer serializer_;
  Trajectory trajectory_;
};

TEST_F(SerializerTest, Test) {
  serializer_.Start(&trajectory_);
  bool done = false;
  do {
    Serializer::Data data = serializer_.Get();
    LOG(ERROR)<<data.size;
    done = data.size == 0;
  } while (!done);
  LOG(ERROR)<<"exiting";
}

}  // namespace base
}  // namespace principia
