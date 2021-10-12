#include "physics/discrete_trajectory_segment_iterator.hpp"

#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/mock_discrete_trajectory_segment.hpp"

namespace principia {
namespace physics {

using base::check_not_null;
using base::not_null;
using geometry::Frame;
using ::testing::Return;

// We use a mock segment in this test to avoid having to go through a
// complicated setup just to test the iterator.
class DiscreteTrajectorySegmentIteratorTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  using Segments = internal_discrete_trajectory_types::Segments<World>;

  DiscreteTrajectorySegmentIterator<World> MakeIterator(
      not_null<Segments const*> const segments,
      Segments::const_iterator const iterator) {
    return DiscreteTrajectorySegmentIterator<World>(segments, iterator);
  }
};

TEST_F(DiscreteTrajectorySegmentIteratorTest, Basic) {
  auto owned_mock1 = std::make_unique<MockDiscreteTrajectorySegment<World>>();
  auto owned_mock2 = std::make_unique<MockDiscreteTrajectorySegment<World>>();
  auto owned_mock3 = std::make_unique<MockDiscreteTrajectorySegment<World>>();
  auto const& mock1 = *owned_mock1;
  auto const& mock2 = *owned_mock2;
  auto const& mock3 = *owned_mock3;

  Segments segments;
  segments.push_back(std::move(*owned_mock1));
  segments.push_back(std::move(*owned_mock2));
  segments.push_back(std::move(*owned_mock3));

  EXPECT_CALL(mock1, size()).WillRepeatedly(Return(5));
  EXPECT_CALL(mock2, size()).WillRepeatedly(Return(1));
  EXPECT_CALL(mock3, size()).WillRepeatedly(Return(3));

  {
    auto iterator = MakeIterator(check_not_null(&segments), segments.begin());
    EXPECT_EQ(5, iterator->size());
    auto const current = ++iterator;
    EXPECT_EQ(1, iterator->size());
    EXPECT_EQ(1, current->size());
    auto const previous = iterator++;
    EXPECT_EQ(3, iterator->size());
    EXPECT_EQ(1, previous->size());
  }
  {
    auto iterator = MakeIterator(check_not_null(&segments), segments.end());
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
