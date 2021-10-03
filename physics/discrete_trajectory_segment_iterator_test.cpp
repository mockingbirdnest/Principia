#include "physics/discrete_trajectory_segment_iterator.hpp"

#include <memory>
#include <vector>

#include "geometry/frame.hpp"
#include "gtest/gtest.h"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/fake_discrete_trajectory.hpp"

namespace principia {
namespace physics {

using geometry::Frame;

class DiscreteTrajectorySegmentIteratorTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  using Segments = internal_discrete_trajectory_types::Segments<World>;

  Segments::value_type NewSegment(std::int64_t size) {
    auto result = std::make_unique<DiscreteTrajectorySegment<World>>();
    for (int i = 0; i < size; ++i) {
      result->timeline.emplace();
    }
    return result;
  }
  DiscreteTrajectorySegmentIterator<World> MakeIterator(
      Segments::const_iterator const iterator) {
    return DiscreteTrajectorySegmentIterator<World>(iterator);
  }
};

TEST_F(DiscreteTrajectorySegmentIteratorTest, Basic) {
  Segments segments;
  segments.push_back(NewSegment(5));
  segments.push_back(NewSegment(1));
  segments.push_back(NewSegment(3));

  {
    auto iterator = MakeIterator(segments.begin());
    EXPECT_EQ(5, iterator->size());
    auto const current = ++iterator;
    EXPECT_EQ(1, iterator->size());
    EXPECT_EQ(1, current->size());
    auto const previous = iterator++;
    EXPECT_EQ(3, iterator->size());
    EXPECT_EQ(1, previous->size());
  }
  {
    auto iterator = MakeIterator(segments.end());
    --iterator;
    EXPECT_EQ(3, (*iterator).size());
    auto const current = --iterator;
    EXPECT_EQ(1, (*iterator).size());
    EXPECT_EQ(1, (*current).size());
    auto const previous = iterator--;
    EXPECT_EQ(5, (*iterator).size());
    EXPECT_EQ(1, (*previous).size());
  }
}

}  // namespace physics
}  // namespace principia
