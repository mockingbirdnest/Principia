#include "physics/discrete_trajectory_segment_iterator.hpp"

#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/mock_discrete_trajectory_segment.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using base::make_not_null_unique;
using base::not_null;
using geometry::Frame;
using geometry::Instant;
using quantities::si::Second;
using ::testing::Return;

// We use a mock segment in this test to avoid having to go through a
// complicated setup just to test the iterator.
class DiscreteTrajectorySegmentIteratorTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;
  using Segments = internal_discrete_trajectory_types::Segments<World>;

  DiscreteTrajectorySegmentIteratorTest() : segments_(MakeSegments(3)) {
    auto it = segments_->begin();
    {
      auto& segment1 = *it;
      segment1.Append(t0_ + 2 * Second, unmoving_origin_);
      segment1.Append(t0_ + 3 * Second, unmoving_origin_);
      segment1.Append(t0_ + 5 * Second, unmoving_origin_);
      segment1.Append(t0_ + 7 * Second, unmoving_origin_);
      segment1.Append(t0_ + 11 * Second, unmoving_origin_);
    }

    ++it;
    {
      auto& segment2 = *it;
      segment2.Append(t0_ + 13 * Second, unmoving_origin_);
    }

    ++it;
    {
      auto& segment3 = *it;
      segment3.Append(t0_ + 13 * Second, unmoving_origin_);
      segment3.Append(t0_ + 17 * Second, unmoving_origin_);
      segment3.Append(t0_ + 19 * Second, unmoving_origin_);
    }
  }

  DiscreteTrajectorySegmentIterator<World> MakeIterator(
      not_null<Segments const*> const segments,
      Segments::const_iterator const iterator) {
    return DiscreteTrajectorySegmentIterator<World>(segments, iterator);
  }

  // Constructs a list of |n| segments which are properly initialized.
  // TODO(phl): Move to a central place.
  static not_null<std::unique_ptr<Segments>> MakeSegments(const int n) {
    auto segments = make_not_null_unique<Segments>(n);
    for (auto it = segments->begin(); it != segments->end(); ++it) {
      *it = DiscreteTrajectorySegment<World>(
          DiscreteTrajectorySegmentIterator<World>(segments.get(), it));
    }
    return segments;
  }

  not_null<std::unique_ptr<Segments>> segments_;
  Instant const t0_;
  DegreesOfFreedom<World> const unmoving_origin_{World::origin,
                                                 World::unmoving};
};

TEST_F(DiscreteTrajectorySegmentIteratorTest, Basic) {
  {
    auto iterator = MakeIterator(segments_.get(), segments_->begin());
    EXPECT_EQ(5, iterator->size());
    auto const current = ++iterator;
    EXPECT_EQ(1, iterator->size());
    EXPECT_EQ(1, current->size());
    auto const previous = iterator++;
    EXPECT_EQ(3, iterator->size());
    EXPECT_EQ(1, previous->size());
  }
  {
    auto iterator = MakeIterator(segments_.get(), segments_->end());
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
