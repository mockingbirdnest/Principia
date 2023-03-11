#include "physics/discrete_trajectory_iterator.hpp"

#include <memory>
#include <optional>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace physics {

using physics::DegreesOfFreedom;
using ::testing::Return;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_named_quantities;
using namespace principia::physics::_discrete_trajectory_iterator;
using namespace principia::physics::_discrete_trajectory_types;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

class DiscreteTrajectoryIteratorTest : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag>;
  using Segments = _discrete_trajectory_types::Segments<World>;

  DiscreteTrajectoryIteratorTest()
      : segments_(MakeSegments(3)) {
    auto it = segments_->begin();
    {
      auto& segment1 = *it;
      EXPECT_OK(segment1.Append(t0_ + 2 * Second, unmoving_origin_));
      EXPECT_OK(segment1.Append(t0_ + 3 * Second, unmoving_origin_));
      EXPECT_OK(segment1.Append(t0_ + 5 * Second, unmoving_origin_));
      EXPECT_OK(segment1.Append(t0_ + 7 * Second, unmoving_origin_));
      EXPECT_OK(segment1.Append(t0_ + 11 * Second, unmoving_origin_));
    }

    ++it;
    {
      auto& segment2 = *it;
      EXPECT_OK(segment2.Append(t0_ + 13 * Second, unmoving_origin_));
    }

    ++it;
    {
      auto& segment3 = *it;
      EXPECT_OK(segment3.Append(t0_ + 13 * Second, unmoving_origin_));
      EXPECT_OK(segment3.Append(t0_ + 17 * Second, unmoving_origin_));
      EXPECT_OK(segment3.Append(t0_ + 19 * Second, unmoving_origin_));
      EXPECT_OK(segment3.Append(t0_ + 23 * Second, unmoving_origin_));
    }
  }

  void Append(Segments::iterator const it,
              Instant const& t,
              DegreesOfFreedom<World> const& degrees_of_freedom) {
    EXPECT_OK(it->Append(t, degrees_of_freedom));
  }

  DiscreteTrajectoryIterator<World> MakeBegin(
      Segments::const_iterator const it) {
    return it->begin();
  }

  DiscreteTrajectoryIterator<World> MakeEnd(
      Segments::const_iterator const it) {
    return it->end();
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

TEST_F(DiscreteTrajectoryIteratorTest, ForwardOneSegment) {
  auto segment = segments_->begin();
  auto iterator = MakeBegin(segment);
  EXPECT_EQ(t0_ + 2 * Second, iterator->time);
  auto const current = ++iterator;
  EXPECT_EQ(t0_ + 3 * Second, iterator->time);
  EXPECT_EQ(t0_ + 3 * Second, current->time);
  auto const previous = iterator++;
  EXPECT_EQ(t0_ + 5 * Second, iterator->time);
  EXPECT_EQ(t0_ + 3 * Second, previous->time);
}

TEST_F(DiscreteTrajectoryIteratorTest, BackwardOneSegment) {
  auto segment = --segments_->end();
  auto iterator = MakeEnd(segment);
  --iterator;
  EXPECT_EQ(t0_ + 23 * Second, (*iterator).time);
  auto const current = --iterator;
  EXPECT_EQ(t0_ + 19 * Second, (*iterator).time);
  EXPECT_EQ(t0_ + 19 * Second, (*current).time);
  auto const previous = iterator--;
  EXPECT_EQ(t0_ + 17 * Second, (*iterator).time);
  EXPECT_EQ(t0_ + 19 * Second, (*previous).time);
}

TEST_F(DiscreteTrajectoryIteratorTest, ForwardAcrossSegments) {
  auto segment = segments_->begin();
  auto iterator = MakeBegin(segment);
  for (int i = 0; i < 4; ++i) {
    ++iterator;
  }
  EXPECT_EQ(t0_ + 11 * Second, iterator->time);
  ++iterator;
  EXPECT_EQ(t0_ + 13 * Second, iterator->time);
  ++iterator;
  EXPECT_EQ(t0_ + 17 * Second, iterator->time);
}

TEST_F(DiscreteTrajectoryIteratorTest, BackwardAcrossSegments) {
  auto segment = --segments_->end();
  auto iterator = MakeEnd(segment);
  for (int i = 0; i < 3; ++i) {
    --iterator;
  }
  EXPECT_EQ(t0_ + 17 * Second, (*iterator).time);
  --iterator;
  EXPECT_EQ(t0_ + 13 * Second, (*iterator).time);
  --iterator;
  EXPECT_EQ(t0_ + 11 * Second, (*iterator).time);
}

TEST_F(DiscreteTrajectoryIteratorTest, Equality) {
  // Construct two iterators that denote the time 13 * Second but in different
  // segments.
  auto it1 = MakeEnd(--segments_->end());
  for (int i = 0; i < 3; ++i) {
    --it1;
  }
  EXPECT_EQ(t0_ + 17 * Second, (*it1).time);
  --it1;
  EXPECT_EQ(t0_ + 13 * Second, (*it1).time);

  auto it2 = MakeBegin(segments_->begin());
  for (int i = 0; i < 4; ++i) {
    ++it2;
  }
  EXPECT_EQ(t0_ + 11 * Second, it2->time);
  ++it2;
  EXPECT_EQ(t0_ + 13 * Second, it2->time);

  EXPECT_EQ(it1, it2);
  EXPECT_NE(it1, MakeBegin(segments_->begin()));
  EXPECT_NE(it2, MakeEnd(--segments_->end()));
  EXPECT_NE(MakeBegin(segments_->begin()), MakeEnd(--segments_->end()));
}

TEST_F(DiscreteTrajectoryIteratorTest, RandomAccess) {
  auto const begin = MakeBegin(segments_->begin());
  auto const end = MakeEnd(--segments_->end());
  {
    auto iterator = begin;
    iterator += 4;
    EXPECT_EQ(t0_ + 11 * Second, iterator->time);
    iterator += 3;
    EXPECT_EQ(t0_ + 19 * Second, iterator->time);
    iterator += -2;
    EXPECT_EQ(t0_ + 13 * Second, iterator->time);
    iterator += 0;
    EXPECT_EQ(t0_ + 13 * Second, iterator->time);
    iterator += 4;
    EXPECT_TRUE(iterator == end);
  }
  {
    auto iterator = begin;
    ++iterator;
    iterator += 0;
    EXPECT_EQ(t0_ + 3 * Second, iterator->time);
    iterator += 6;
    EXPECT_EQ(t0_ + 19 * Second, iterator->time);
    iterator += -5;
    EXPECT_EQ(t0_ + 5 * Second, iterator->time);
    iterator += 3;
    EXPECT_EQ(t0_ + 13 * Second, iterator->time);
    iterator += 4;
    EXPECT_TRUE(iterator == end);
  }
  {
    auto iterator = MakeEnd(--segments_->end());
    iterator -= 4;
    EXPECT_EQ(t0_ + 13 * Second, iterator->time);
    iterator -= -3;
    EXPECT_EQ(t0_ + 23 * Second, iterator->time);
    iterator -= 2;
    EXPECT_EQ(t0_ + 17 * Second, iterator->time);
    iterator -= 0;
    EXPECT_EQ(t0_ + 17 * Second, iterator->time);
    iterator -= 6;
    EXPECT_TRUE(iterator == begin);
  }
  {
    auto iterator = MakeEnd(--segments_->end());
    --iterator;
    --iterator;
    iterator -= 0;
    EXPECT_EQ(t0_ + 19 * Second, iterator->time);
    iterator -= 6;
    EXPECT_EQ(t0_ + 3 * Second, iterator->time);
    iterator -= -5;
    EXPECT_EQ(t0_ + 17 * Second, iterator->time);
    iterator -= 1;
    EXPECT_EQ(t0_ + 13 * Second, iterator->time);
    iterator -= 5;
    EXPECT_TRUE(iterator == begin);
  }
  {
    auto iterator = begin;
    ++iterator;
    ++iterator;
    EXPECT_EQ(t0_ + 3 * Second, (iterator - 1)->time);
    EXPECT_EQ(t0_ + 19 * Second, (iterator + 5)->time);
    EXPECT_EQ(t0_ + 13 * Second, (3 + iterator)->time);
    EXPECT_EQ(t0_ + 7 * Second, iterator[1].time);
    EXPECT_EQ(t0_ + 2 * Second, iterator[-2].time);
  }
  {
    auto iterator1 = begin;
    ++iterator1;
    auto const iterator2 = iterator1 + 6;
    EXPECT_EQ(6, iterator2 - iterator1);
    EXPECT_EQ(0, iterator2 - iterator2);
    auto const iterator3 = end;
    EXPECT_EQ(2, iterator3 - iterator2);
    EXPECT_EQ(8, iterator3 - iterator1);
  }
  {
    EXPECT_LT(begin + 3, begin + 4);
    EXPECT_LT(begin + 4, end);
    EXPECT_LE(begin + 4, begin + 4);
    EXPECT_LE(begin + 3, begin + 4);
    EXPECT_LE(begin + 4, end);
    EXPECT_LE(end, end);
    EXPECT_GE(begin + 4, begin + 4);
    EXPECT_GE(begin + 4, begin + 3);
    EXPECT_GE(end, begin + 4);
    EXPECT_GE(end, end);
    EXPECT_GT(begin + 4, begin + 3);
    EXPECT_GT(end, begin + 4);
  }
}

// Empty segments may exist in a transient manner, we must be able to iterate
// over them.
TEST_F(DiscreteTrajectoryIteratorTest, EmptySegment) {
  auto segments = MakeSegments(1);
  {
    int count = 0;
    for (auto const& point : segments->front()) {
      ++count;
    }
    EXPECT_EQ(0, count);
  }
  {
    int count = 0;
    for (auto it = segments->front().rbegin();
         it != segments->front().rend();
         ++it) {
      ++count;
    }
    EXPECT_EQ(0, count);
  }
}

// Check that repeated points don't cause confusion regarding the end of a
// segment.
TEST_F(DiscreteTrajectoryIteratorTest, SegmentEnd) {
  auto segment0 = segments_->begin();
  auto segment1 = std::next(segment0);
  auto iterator = segment0->begin();
  for (int i = 0; i < 5; ++i) {
    ++iterator;
  }
  EXPECT_TRUE(iterator != segment1->end());
}

// Checkt that rbegin() works if the next segment is empty.
TEST_F(DiscreteTrajectoryIteratorTest, EmptyLastSegment) {
  auto segments = MakeSegments(2);
  auto segment = segments->begin();
  Append(segment, t0_, unmoving_origin_);
  EXPECT_EQ(t0_, segment->rbegin()->time);
}

}  // namespace physics
}  // namespace principia
