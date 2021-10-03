#include "physics/discrete_trajectory_segment_iterator.hpp"

#include <memory>
#include <vector>

#include "geometry/frame.hpp"
#include "gtest/gtest.h"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {

using geometry::Frame;

namespace {

// Fake the segment class to avoid complicated dependencies.
template<typename Frame>
class FakeDiscreteTrajectorySegment : public DiscreteTrajectorySegment<Frame> {
 public:
  FakeDiscreteTrajectorySegment(std::vector<int> const& segment);

  virtual std::int64_t size() const override;

 private:
  std::vector<int> segment_;
};

template<typename Frame>
FakeDiscreteTrajectorySegment<Frame>::FakeDiscreteTrajectorySegment(
    std::vector<int> const& segment)
    : segment_(segment) {}

template<typename Frame>
std::int64_t FakeDiscreteTrajectorySegment<Frame>::size() const {
  return segment_.size();
}

}  // namespace

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
  Segments segments;
  segments.push_back(
      std::make_unique<FakeDiscreteTrajectorySegment<World>>(
          std::vector{2, 3, 5, 7, 11}));
  segments.push_back(
      std::make_unique<FakeDiscreteTrajectorySegment<World>>(
          std::vector{13}));
  segments.push_back(
      std::make_unique<FakeDiscreteTrajectorySegment<World>>(
          std::vector{17, 19, 23}));

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
