#include "physics/discrete_trajectory_segment_iterator.hpp"

#include <memory>
#include <vector>

#include "gtest/gtest.h"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {

// Fake the segment class to avoid complicated dependencies.
template<typename Frame>
class DiscreteTrajectorySegment {
 public:
  std::vector<int> segment;
};

namespace internal_discrete_trajectory_segment_iterator {

class DiscreteTrajectorySegmentIteratorTest : public ::testing::Test {
 protected:
  using Frame = void;

  using Segments = internal_discrete_trajectory_types::Segments<Frame>;

  DiscreteTrajectorySegmentIterator<Frame> MakeIterator(
    Segments::const_iterator const iterator) {
    return DiscreteTrajectorySegmentIterator<Frame>(iterator);
  }
};

TEST_F(DiscreteTrajectorySegmentIteratorTest, Basic) {
  DiscreteTrajectorySegment<Frame> const primes1{{2, 3, 5, 7, 11}};
  DiscreteTrajectorySegment<Frame> const primes2{{13}};
  DiscreteTrajectorySegment<Frame> const primes3{{17, 19, 23}};

  Segments segments;
  segments.push_back(
      std::make_unique<DiscreteTrajectorySegment<Frame>>(primes1));
  segments.push_back(
      std::make_unique<DiscreteTrajectorySegment<Frame>>(primes2));
  segments.push_back(
      std::make_unique<DiscreteTrajectorySegment<Frame>>(primes3));

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

}  // namespace internal_discrete_trajectory_segment_iterator
}  // namespace physics
}  // namespace principia
