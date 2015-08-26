#include "physics/forkable.hpp"

#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"

namespace principia {

using geometry::Instant;

namespace physics {

class FakeTrajectory : public Forkable<FakeTrajectory> {
 public:
  using TimelineConstIterator = std::vector<Instant>::const_iterator;

 protected:
  TimelineConstIterator timeline_end() const override;
  TimelineConstIterator timeline_find(Instant const& time) const override;
  void timeline_insert(TimelineConstIterator begin,
                       TimelineConstIterator end) override;
  bool timeline_empty() const override;

 private:
  std::vector<Instant> timeline_;
};

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
}

void FakeTrajectory::timeline_insert(TimelineConstIterator begin,
                                     TimelineConstIterator end) {
  CHECK(timeline_empty());
  timeline_.insert(begin, end);
}

bool FakeTrajectory::timeline_empty() const {
  return timeline_.empty();
}

class ForkableTest : public testing::Test {
};

}  // namespace physics
}  // namespace principia
