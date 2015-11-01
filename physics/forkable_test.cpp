#include "physics/forkable.hpp"

#include <list>
#include <vector>

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {

using geometry::Instant;
using quantities::si::Second;
using ::testing::ElementsAre;

namespace physics {

class FakeTrajectory;

template<>
struct ForkableTraits<FakeTrajectory> {
  using TimelineConstIterator = std::list<Instant>::const_iterator;
  static Instant const& time(TimelineConstIterator const it);
};

class FakeTrajectory : public Forkable<FakeTrajectory> {
 public:
  FakeTrajectory() = default;

  void push_back(Instant const& time);

  using Forkable<FakeTrajectory>::NewFork;
  using Forkable<FakeTrajectory>::DeleteAllForksAfter;
  using Forkable<FakeTrajectory>::DeleteAllForksBefore;

 protected:
  not_null<FakeTrajectory*> that() override;
  not_null<FakeTrajectory const*> that() const override;

  TimelineConstIterator timeline_begin() const override;
  TimelineConstIterator timeline_end() const override;
  TimelineConstIterator timeline_find(Instant const& time) const override;
  TimelineConstIterator timeline_lower_bound(
                            Instant const& time) const override;
  bool timeline_empty() const override;

 private:
  // Use list<> because we want the iterators to remain valid across operations.
  std::list<Instant> timeline_;

  template<typename Tr4jectory>
  friend class Forkable;
};

Instant const& ForkableTraits<FakeTrajectory>::time(
  TimelineConstIterator const it) {
  return *it;
}

void FakeTrajectory::push_back(Instant const& time) {
  timeline_.push_back(time);
}

not_null<FakeTrajectory*> FakeTrajectory::that() {
  return this;
}

not_null<FakeTrajectory const*> FakeTrajectory::that() const {
  return this;
}

FakeTrajectory::TimelineConstIterator FakeTrajectory::timeline_begin() const {
  return timeline_.begin();
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

FakeTrajectory::TimelineConstIterator FakeTrajectory::timeline_lower_bound(
                                          Instant const& time) const {
  // Stupid O(N) search.
  for (auto it = timeline_.begin(); it != timeline_.end(); ++it) {
    if (*it >= time) {
      return it;
    }
  }
  return timeline_.end();
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

  static std::vector<Instant> After(
      not_null<FakeTrajectory const*> const trajectory,
      Instant const& time) {
    std::vector<Instant> after;
    for (FakeTrajectory::Iterator it = trajectory->Find(time);
         it != trajectory->End();
         ++it) {
      after.push_back(*it.current());
    }
    return after;
  }

  static Instant const& LastTime(
      not_null<FakeTrajectory const*> const trajectory) {
    FakeTrajectory::Iterator it = trajectory->End();
    --it;
    return *it.current();
  }

  static std::vector<Instant> Times(
      not_null<FakeTrajectory const*> const trajectory) {
    std::vector<Instant> times;
    for (FakeTrajectory::Iterator it = trajectory->Begin();
         it != trajectory->End();
         ++it) {
      times.push_back(*it.current());
    }
    return times;
  }

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

TEST_F(ForkableTest, ForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
  fork->push_back(t4_);
  auto times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  times = Times(fork);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t4_));
}

TEST_F(ForkableTest, ForkAtLast) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork1 = trajectory_.NewFork(t3_);
  not_null<FakeTrajectory*> const fork2 = fork1->NewFork(LastTime(fork1));
  not_null<FakeTrajectory*> const fork3 = fork2->NewFork(LastTime(fork1));
  EXPECT_EQ(t3_, LastTime(&trajectory_));
  EXPECT_EQ(t3_, LastTime(fork1));

  auto times = Times(fork2);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(t3_, LastTime(fork2));
  EXPECT_EQ(t3_, *fork2->Fork().current());

  auto after = After(fork3, t3_);
  EXPECT_THAT(after, ElementsAre(t3_));

  after = After(fork2, t3_);
  EXPECT_THAT(after, ElementsAre(t3_));

  fork1->push_back(t4_);
  times = Times(fork2);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  after = After(fork1, t3_);
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  times = Times(fork3);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  fork2->push_back(t4_);
  after = After(fork2, t3_);
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  fork3->push_back(t4_);
  after = After(fork3, t3_);
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  after = After(fork3, t2_);
  EXPECT_THAT(after, ElementsAre(t2_, t3_, t4_));
}

TEST_F(ForkableDeathTest, DeleteForkError) {
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    FakeTrajectory* root = &trajectory_;
    trajectory_.DeleteFork(&root);
  }, "fork.*End");
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    FakeTrajectory* fork1 = trajectory_.NewFork(t1_);
    fork1->push_back(t2_);
    FakeTrajectory* fork2 = fork1->NewFork(t2_);
    trajectory_.DeleteFork(&fork2);
  }, "not a child");
}

TEST_F(ForkableTest, DeleteForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork1 = trajectory_.NewFork(t2_);
  FakeTrajectory* fork2 = trajectory_.NewFork(t2_);
  fork1->push_back(t4_);
  trajectory_.DeleteFork(&fork2);
  EXPECT_EQ(nullptr, fork2);
  auto times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  times = Times(fork1);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t4_));
}

TEST_F(ForkableDeathTest, DeleteAllForksAfterError) {
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    trajectory_.push_back(t2_);
    not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
    fork->DeleteAllForksAfter(t1_);
  }, "before the fork time");
}

TEST_F(ForkableTest, DeleteAllForksAfterSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
  fork->push_back(t4_);

  fork->DeleteAllForksAfter(t3_ + (t4_ - t3_) / 2);
  auto times = Times(fork);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t4_));

  fork->DeleteAllForksAfter(t2_);
  times = Times(fork);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t4_));

  times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  trajectory_.DeleteAllForksAfter(t1_);
  times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  // Don't use fork, it is dangling.
}

TEST_F(ForkableDeathTest, DeleteAllForksBeforeError) {
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t1_);
    fork->DeleteAllForksBefore(t1_);
  }, "nonroot");
}

TEST_F(ForkableTest, DeleteAllForksBeforeSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
  fork->push_back(t4_);

  trajectory_.DeleteAllForksBefore(t1_ + (t2_ - t1_) / 2);
  auto times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  times = Times(fork);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t4_));

  trajectory_.DeleteAllForksBefore(t2_);
  times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  // Don't use fork, it is dangling.
}

TEST_F(ForkableDeathTest, IteratorDecrementError) {
  EXPECT_DEATH({
    auto it = trajectory_.End();
    --it;
  }, "parent_.*non NULL");
}

TEST_F(ForkableTest, IteratorDecrementNoForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  auto it = trajectory_.End();
  --it;
  EXPECT_EQ(t3_, *it.current());
  --it;
  EXPECT_EQ(t2_, *it.current());
  --it;
  EXPECT_EQ(t1_, *it.current());
}

TEST_F(ForkableTest, IteratorDecrementForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  auto fork = trajectory_.NewFork(t1_);
  trajectory_.push_back(t4_);
  fork->push_back(t3_);
  auto it = fork->End();
  --it;
  EXPECT_EQ(t3_, *it.current());
  --it;
  EXPECT_EQ(t1_, *it.current());
}

TEST_F(ForkableTest, IteratorDecrementMultipleForksSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  auto fork1 = trajectory_.NewFork(t2_);
  auto fork2 = fork1->NewFork(t2_);
  auto fork3 = fork2->NewFork(t2_);
  fork2->push_back(t3_);
  auto it = fork3->End();
  --it;
  EXPECT_EQ(t2_, *it.current());
  --it;
  EXPECT_EQ(t1_, *it.current());
  EXPECT_EQ(it, fork3->Begin());
}

TEST_F(ForkableDeathTest, IteratorIncrementError) {
  EXPECT_DEATH({
    auto it = trajectory_.Begin();
    ++it;
  }, "current.*!=.*end");
}

TEST_F(ForkableTest, IteratorIncrementNoForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  auto it = trajectory_.Begin();
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(t3_, *it.current());
}

TEST_F(ForkableTest, IteratorIncrementForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  auto fork = trajectory_.NewFork(t1_);
  trajectory_.push_back(t4_);
  fork->push_back(t3_);
  auto it = fork->Begin();
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t3_, *it.current());
  ++it;
  EXPECT_EQ(it, fork->End());
}

TEST_F(ForkableTest, IteratorIncrementMultipleForksSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  auto fork1 = trajectory_.NewFork(t2_);
  auto fork2 = fork1->NewFork(t2_);
  auto fork3 = fork2->NewFork(t2_);
  auto it = fork3->Begin();
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(it, fork3->End());
  fork3->push_back(t3_);
  --it;
  EXPECT_EQ(t3_, *it.current());
  it = fork3->Begin();
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(t3_, *it.current());
  ++it;
  EXPECT_EQ(it, fork3->End());
}

#if !defined(_DEBUG)
TEST_F(ForkableTest, IteratorEndEquality) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  auto fork1 = trajectory_.NewFork(t1_);
  auto fork2 = trajectory_.NewFork(t2_);
  auto it1 = fork1->End();
  auto it2 = fork2->End();
  EXPECT_NE(it1, it2);
}
#endif

TEST_F(ForkableTest, Root) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
  EXPECT_TRUE(trajectory_.is_root());
  EXPECT_FALSE(fork->is_root());
  EXPECT_EQ(&trajectory_, trajectory_.root());
  EXPECT_EQ(&trajectory_, fork->root());
  EXPECT_TRUE(trajectory_.Fork() == trajectory_.End());
  EXPECT_EQ(t2_, *fork->Fork().current());
}

TEST_F(ForkableTest, IteratorBeginSuccess) {
  auto it = trajectory_.Begin();
  EXPECT_EQ(it, trajectory_.End());

  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);

  it = trajectory_.Begin();
  EXPECT_NE(it, trajectory_.End());
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(t3_, *it.current());
  ++it;
  EXPECT_EQ(it, trajectory_.End());

  not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
  fork->push_back(t4_);

  it = fork->Begin();
  EXPECT_NE(it, fork->End());
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(t4_, *it.current());
  ++it;
  EXPECT_EQ(it, fork->End());
}

TEST_F(ForkableTest, IteratorFindSuccess) {
  auto it = trajectory_.Find(t0_);
  EXPECT_EQ(it, trajectory_.End());

  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);

  it = trajectory_.Find(t0_);
  EXPECT_EQ(it, trajectory_.End());
  it = trajectory_.Find(t1_);
  EXPECT_NE(it, trajectory_.End());
  EXPECT_EQ(t1_, *it.current());
  it = trajectory_.Find(t2_);
  EXPECT_EQ(t2_, *it.current());
  it = trajectory_.Find(t4_);
  EXPECT_EQ(it, trajectory_.End());

  not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
  fork->push_back(t4_);

  it = fork->Find(t0_);
  EXPECT_EQ(it, fork->End());
  it = fork->Find(t1_);
  EXPECT_NE(it, fork->End());
  EXPECT_EQ(t1_, *it.current());
  it = fork->Find(t2_);
  EXPECT_EQ(t2_, *it.current());
  it = fork->Find(t4_);
  EXPECT_EQ(t4_, *it.current());
  it = fork->Find(t4_ + 1 * Second);
  EXPECT_EQ(it, fork->End());
}

TEST_F(ForkableTest, IteratorLowerBoundSuccess) {
  auto it = trajectory_.LowerBound(t0_);
  EXPECT_EQ(it, trajectory_.End());

  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);

  it = trajectory_.LowerBound(t0_);
  EXPECT_EQ(t1_, *it.current());
  it = trajectory_.LowerBound(t1_);
  EXPECT_EQ(t1_, *it.current());
  it = trajectory_.LowerBound(t2_);
  EXPECT_EQ(t2_, *it.current());
  it = trajectory_.LowerBound(t4_);
  EXPECT_EQ(it, trajectory_.End());

  not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
  fork->push_back(t4_);

  it = fork->LowerBound(t0_);
  EXPECT_EQ(t1_, *it.current());
  it = fork->LowerBound(t1_);
  EXPECT_NE(it, fork->End());
  EXPECT_EQ(t1_, *it.current());
  it = fork->LowerBound(t2_);
  EXPECT_EQ(t2_, *it.current());
  it = fork->LowerBound(t4_);
  EXPECT_EQ(t4_, *it.current());
  it = fork->LowerBound(t4_ + 1 * Second);
  EXPECT_EQ(it, fork->End());
}

TEST_F(ForkableTest, IteratorSerializationSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork1 = trajectory_.NewFork(t2_);
  not_null<FakeTrajectory*> const fork2 = trajectory_.NewFork(t2_);
  not_null<FakeTrajectory*> const fork3 = trajectory_.NewFork(t3_);
  fork2->DeleteAllForksAfter(t2_);
  fork3->push_back(t4_);

  {
    serialization::Trajectory::Pointer message;
    trajectory_.WritePointerToMessage(&message);
    EXPECT_EQ(0, message.fork_size());
    auto trajectory = FakeTrajectory::ReadPointerFromMessage(message,
                                                             &trajectory_);
    EXPECT_EQ(&trajectory_, trajectory);
    auto it = trajectory->Begin();
    EXPECT_EQ(t1_, *it.current());
    ++it;
    EXPECT_EQ(t2_, *it.current());
    ++it;
    EXPECT_EQ(t3_, *it.current());
  }

  {
    serialization::Trajectory::Pointer message;
    fork1->WritePointerToMessage(&message);
    EXPECT_EQ(1, message.fork_size());
    EXPECT_EQ(0, message.fork(0).children_distance());
    EXPECT_EQ(1, message.fork(0).timeline_distance());
    auto trajectory = FakeTrajectory::ReadPointerFromMessage(message,
                                                             &trajectory_);
    EXPECT_EQ(fork1, trajectory);
    auto it = trajectory->Begin();
    EXPECT_EQ(t1_, *it.current());
    ++it;
    EXPECT_EQ(t2_, *it.current());
    ++it;
    EXPECT_EQ(it, trajectory->End());
  }

  {
    serialization::Trajectory::Pointer message;
    fork2->WritePointerToMessage(&message);
    EXPECT_EQ(1, message.fork_size());
    EXPECT_EQ(1, message.fork(0).children_distance());
    EXPECT_EQ(1, message.fork(0).timeline_distance());
    auto trajectory = FakeTrajectory::ReadPointerFromMessage(message,
                                                             &trajectory_);
    EXPECT_EQ(fork2, trajectory);
    auto it = trajectory->Begin();
    EXPECT_EQ(t1_, *it.current());
    ++it;
    EXPECT_EQ(t2_, *it.current());
    ++it;
    EXPECT_EQ(it, trajectory->End());
  }

  {
    serialization::Trajectory::Pointer message;
    fork3->WritePointerToMessage(&message);
    EXPECT_EQ(1, message.fork_size());
    EXPECT_EQ(2, message.fork(0).children_distance());
    EXPECT_EQ(2, message.fork(0).timeline_distance());
    auto trajectory = FakeTrajectory::ReadPointerFromMessage(message,
                                                             &trajectory_);
    EXPECT_EQ(fork3, trajectory);
    auto it = trajectory->Begin();
    EXPECT_EQ(t1_, *it.current());
    ++it;
    EXPECT_EQ(t2_, *it.current());
    ++it;
    EXPECT_EQ(t3_, *it.current());
    ++it;
    EXPECT_EQ(t4_, *it.current());
    ++it;
    EXPECT_EQ(it, trajectory->End());
  }
}

}  // namespace physics
}  // namespace principia
