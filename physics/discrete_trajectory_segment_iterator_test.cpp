#include "physics/discrete_trajectory_segment_iterator.hpp"

#include <memory>
#include <vector>

#include "geometry/frame.hpp"
#include "gtest/gtest.h"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {

using geometry::Frame;

// Fake the segment class to avoid complicated dependencies.
namespace {
template<typename Frame>
class DiscreteTrajectorySegment {
 public:
  std::vector<int> segment;
};
}

class DiscreteTrajectorySegmentIteratorTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  using Segments = internal_discrete_trajectory_types::Segments<World>;

  DiscreteTrajectorySegmentIterator<World> MakeIterator(
      Segments::const_iterator const iterator) {
    return DiscreteTrajectorySegmentIterator<World>(iterator);
  }
};

TEST_F(DiscreteTrajectorySegmentIteratorTest, Basic) {
  DiscreteTrajectorySegment<World> const primes1{{2, 3, 5, 7, 11}};
  DiscreteTrajectorySegment<World> const primes2{{13}};
  DiscreteTrajectorySegment<World> const primes3{{17, 19, 23}};

  Segments segments;
  segments.push_back(
      std::make_unique<DiscreteTrajectorySegment<World>>(primes1));
  segments.push_back(
      std::make_unique<DiscreteTrajectorySegment<World>>(primes2));
  segments.push_back(
      std::make_unique<DiscreteTrajectorySegment<World>>(primes3));

  {
    auto iterator = MakeIterator(segments.begin());
    EXPECT_EQ(5, iterator->segment.size());
    auto const current = ++iterator;
    EXPECT_EQ(1, iterator->segment.size());
    EXPECT_EQ(1, current->segment.size());
    auto const previous = iterator++;
    EXPECT_EQ(3, iterator->segment.size());
    EXPECT_EQ(1, previous->segment.size());
  }
  {
    auto iterator = MakeIterator(segments.end());
    --iterator;
    EXPECT_EQ(3, (*iterator).segment.size());
    auto const current = --iterator;
    EXPECT_EQ(1, (*iterator).segment.size());
    EXPECT_EQ(1, (*current).segment.size());
    auto const previous = iterator--;
    EXPECT_EQ(5, (*iterator).segment.size());
    EXPECT_EQ(1, (*previous).segment.size());
  }
}

}  // namespace physics
}  // namespace principia
