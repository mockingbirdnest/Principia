#include "physics/forkable.hpp"

#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {

using geometry::Instant;
using si::Second;

namespace physics {

class FakeTrajectory : public Forkable<std::vector<Instant>::const_iterator> {
 public:
  FakeTrajectory();

  void push_back(Instant const& time);

 protected:
  TimelineConstIterator timeline_end() const override;
  TimelineConstIterator timeline_find(Instant const& time) const override;
  void timeline_insert(TimelineConstIterator begin,
                       TimelineConstIterator end) override;
  bool timeline_empty() const override;

 private:
  std::vector<Instant> timeline_;
};

FakeTrajectory::FakeTrajectory()
    : Forkable<std::vector<Instant>::const_iterator>() {}

void FakeTrajectory::push_back(Instant const& time) {
  timeline_.push_back(time);
}

FakeTrajectory::TimelineConstIterator FakeTrajectory::timeline_end() const {
  return timeline_.end();
}

FakeTrajectory::TimelineConstIterator FakeTrajectory::timeline_find(
    Instant const & time) const {
  // Stupid O(N) search.
  for (auto it = timeline_.begin(); it != timeline_.end(); ++it) {
    if (*it == time) {
      return it;
    }
  }
  return timeline_.end();
}

void FakeTrajectory::timeline_insert(TimelineConstIterator begin,
                                     TimelineConstIterator end) {
  CHECK(timeline_empty());
  timeline_.insert(timeline_.end(), begin, end);
}

bool FakeTrajectory::timeline_empty() const {
  return timeline_.empty();
}

class ForkableTest : public testing::Test {
 protected:
  ForkableTest() :
    t0_(),
    t1_(t0_ + 7 * Second),
    t2_(t0_ + 17 * Second),
    t3_(t0_ + 27 * Second),
    t4_(t0_ + 37 * Second) {}

  FakeTrajectory trajectory_;
  Instant t0_, t1_, t2_, t3_, t4_;
};

using ForkableDeathTest = ForkableTest;

TEST_F(ForkableDeathTest, ForkError) {
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    trajectory_.push_back(t3_);
    trajectory_.NewFork(t2_);
  }, "nonexistent time");
}

}  // namespace physics
}  // namespace principia
