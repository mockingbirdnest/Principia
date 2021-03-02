﻿
#include "physics/forkable.hpp"

#include <list>
#include <vector>

#include "base/not_constructible.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_forkable {

using base::not_constructible;
using geometry::Instant;
using quantities::si::Second;
using ::testing::ElementsAre;

struct FakeTrajectoryTraits : not_constructible {
  using TimelineConstIterator = std::list<Instant>::const_iterator;
  static Instant const& time(TimelineConstIterator const it);
};

class FakeTrajectory;

class FakeTrajectoryIterator : public ForkableIterator<FakeTrajectory,
                                                       FakeTrajectoryIterator,
                                                       FakeTrajectoryTraits> {
 public:
  using ForkableIterator<FakeTrajectory,
                         FakeTrajectoryIterator,
                         FakeTrajectoryTraits>::current;
  using reference = Instant const&;

  reference operator*() const;

 protected:
  not_null<FakeTrajectoryIterator*> that() override;
  not_null<FakeTrajectoryIterator const*> that() const override;
};

class FakeTrajectory : public Forkable<FakeTrajectory,
                                       FakeTrajectoryIterator,
                                       FakeTrajectoryTraits> {
 public:
  using Iterator = FakeTrajectoryIterator;

  FakeTrajectory() = default;

  void pop_front();
  void push_front(Instant const& time);
  void push_back(Instant const& time);

  using Forkable<FakeTrajectory, Iterator, FakeTrajectoryTraits>::
      NewFork;
  using Forkable<FakeTrajectory, Iterator, FakeTrajectoryTraits>::
      AttachForkToCopiedBegin;
  using Forkable<FakeTrajectory, Iterator, FakeTrajectoryTraits>::
      DetachForkWithCopiedBegin;
  using Forkable<FakeTrajectory, Iterator, FakeTrajectoryTraits>::
      DeleteAllForksAfter;
  using Forkable<FakeTrajectory, Iterator, FakeTrajectoryTraits>::
      CheckNoForksBefore;

  TimelineConstIterator timeline_begin() const override;
  TimelineConstIterator timeline_end() const override;
  TimelineConstIterator timeline_find(Instant const& time) const override;
  TimelineConstIterator timeline_lower_bound(
                            Instant const& time) const override;
  bool timeline_empty() const override;
  std::int64_t timeline_size() const override;

 protected:
  not_null<FakeTrajectory*> that() override;
  not_null<FakeTrajectory const*> that() const override;

 private:
  // Use list<> because we want the iterators to remain valid across operations.
  std::list<Instant> timeline_;

  template<typename, typename, typename>
  friend class ForkableIterator;
  template<typename, typename, typename>
  friend class Forkable;
};

Instant const& FakeTrajectoryTraits::time(TimelineConstIterator const it) {
  return *it;
}

FakeTrajectoryIterator::reference FakeTrajectoryIterator::operator*() const {
  return *current();
}

not_null<FakeTrajectoryIterator*> FakeTrajectoryIterator::that() {
  return this;
}

not_null<FakeTrajectoryIterator const*> FakeTrajectoryIterator::that() const {
  return this;
}

void FakeTrajectory::pop_front() {
  timeline_.pop_front();
}

void FakeTrajectory::push_front(Instant const& time) {
  timeline_.push_front(time);
}

void FakeTrajectory::push_back(Instant const& time) {
  timeline_.push_back(time);
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

std::int64_t FakeTrajectory::timeline_size() const {
  return timeline_.size();
}

not_null<FakeTrajectory*> FakeTrajectory::that() {
  return this;
}

not_null<FakeTrajectory const*> FakeTrajectory::that() const {
  return this;
}

class ForkableTest : public testing::Test {
 protected:
  ForkableTest() :
    t0_(),
    t1_(t0_ + 7 * Second),
    t2_(t0_ + 17 * Second),
    t3_(t0_ + 27 * Second),
    t4_(t0_ + 37 * Second),
    t5_(t0_ + 47 * Second) {}

  static std::vector<Instant> After(
      not_null<FakeTrajectory const*> const trajectory,
      Instant const& time) {
    std::vector<Instant> after;
    for (FakeTrajectory::Iterator it = trajectory->Find(time);
         it != trajectory->end();
         ++it) {
      after.push_back(*it.current());
    }
    return after;
  }

  static Instant const& LastTime(
      not_null<FakeTrajectory const*> const trajectory) {
    FakeTrajectory::Iterator it = trajectory->end();
    --it;
    return *it.current();
  }

  static std::vector<Instant> Times(
      not_null<FakeTrajectory const*> const trajectory) {
    std::vector<Instant> times;
    for (FakeTrajectory::Iterator it = trajectory->begin();
         it != trajectory->end();
         ++it) {
      times.push_back(*it.current());
    }
    return times;
  }

  FakeTrajectory trajectory_;
  Instant t0_, t1_, t2_, t3_, t4_, t5_;
};

using ForkableDeathTest = ForkableTest;

TEST_F(ForkableDeathTest, ForkError) {
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    trajectory_.push_back(t3_);
    trajectory_.NewFork(trajectory_.timeline_find(t2_));
  }, "!is_root");
  EXPECT_DEATH({
    trajectory_.Fork();
  }, "!is_root");
}

TEST_F(ForkableTest, ForkSuccess) {
  EXPECT_TRUE(trajectory_.Empty());
  trajectory_.push_back(t1_);
  EXPECT_EQ(1, trajectory_.Size());
  EXPECT_FALSE(trajectory_.Empty());
  trajectory_.push_back(t2_);
  EXPECT_EQ(2, trajectory_.Size());
  EXPECT_FALSE(trajectory_.Empty());
  trajectory_.push_back(t3_);
  EXPECT_EQ(3, trajectory_.Size());
  EXPECT_FALSE(trajectory_.Empty());
  not_null<FakeTrajectory*> const fork =
       trajectory_.NewFork(trajectory_.timeline_find(t2_));
  EXPECT_EQ(2, fork->Size());
  EXPECT_FALSE(fork->Empty());
  fork->push_back(t4_);
  EXPECT_EQ(3, fork->Size());
  EXPECT_FALSE(fork->Empty());
  auto times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  times = Times(fork);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t4_));
}

TEST_F(ForkableTest, Size) {
  EXPECT_TRUE(trajectory_.Empty());
  trajectory_.push_back(t1_);
  EXPECT_EQ(1, trajectory_.Size());
  not_null<FakeTrajectory*> const fork1 =
      trajectory_.NewFork(trajectory_.timeline_find(t1_));
  EXPECT_EQ(1, fork1->Size());
  not_null<FakeTrajectory*> const fork2 = fork1->NewFork(fork1->timeline_end());
  fork2->push_back(t2_);
  EXPECT_EQ(2, fork2->Size());
}

TEST_F(ForkableTest, ForkAtLast) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork1 =
      trajectory_.NewFork(trajectory_.timeline_find(t3_));
  not_null<FakeTrajectory*> const fork2 =
      fork1->NewFork(fork1->timeline_find(LastTime(fork1)));
  not_null<FakeTrajectory*> const fork3 =
      fork2->NewFork(fork2->timeline_find(LastTime(fork1)));
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
    trajectory_.DeleteFork(root);
  }, "!is_root");
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    FakeTrajectory* fork1 = trajectory_.NewFork(trajectory_.timeline_find(t1_));
    fork1->push_back(t2_);
    FakeTrajectory* fork2 = fork1->NewFork(fork1->timeline_find(t2_));
    trajectory_.DeleteFork(fork2);
  }, "not a child");
}

TEST_F(ForkableTest, DeleteForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork1 =
      trajectory_.NewFork(trajectory_.timeline_find(t2_));
  FakeTrajectory* fork2 = trajectory_.NewFork(trajectory_.timeline_find(t2_));
  fork1->push_back(t4_);
  trajectory_.DeleteFork(fork2);
  EXPECT_EQ(nullptr, fork2);
  auto times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  times = Times(fork1);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t4_));
}

TEST_F(ForkableDeathTest, AttachForkWithCopiedBeginError) {
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    Box<FakeTrajectory> fork{};
    trajectory_.AttachForkToCopiedBegin(std::move(fork));
  }, "timeline_empty");
}

TEST_F(ForkableTest, AttachForkWithCopiedBeginSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);

  Box<FakeTrajectory> fork1{};
  fork1->push_back(t3_);
  not_null<FakeTrajectory*> const fork2 =
      fork1->NewFork(fork1->timeline_find(t3_));
  fork2->push_back(t4_);
  auto times = Times(fork1.get());
  EXPECT_THAT(times, ElementsAre(t3_));
  times = Times(fork2);
  EXPECT_THAT(times, ElementsAre(t3_, t4_));

  auto const unowned_fork1 = fork1.get();
  trajectory_.AttachForkToCopiedBegin(std::move(fork1));
  unowned_fork1->pop_front();

  times = Times(unowned_fork1);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  times = Times(fork2);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_, t4_));
}

TEST_F(ForkableTest, AttachForkWithCopiedBeginEmpty) {
  trajectory_.push_back(t1_);
  not_null<FakeTrajectory*> const fork1 =
      trajectory_.NewFork(trajectory_.timeline_find(t1_));
  Box<FakeTrajectory> fork2{};
  fork2->push_back(t3_);
  fork1->AttachForkToCopiedBegin(std::move(fork2));
}

TEST_F(ForkableDeathTest, DetachForkWithCopiedBeginError) {
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    trajectory_.DetachForkWithCopiedBegin();
  }, "!is_root");
}

TEST_F(ForkableTest, DetachForkWithCopiedBeginSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork1 =
      trajectory_.NewFork(trajectory_.timeline_find(t2_));
  FakeTrajectory* fork2 = trajectory_.NewFork(trajectory_.timeline_find(t2_));
  FakeTrajectory* fork3 = fork1->NewFork(fork1->timeline_find(t2_));
  fork1->push_back(t4_);

  fork1->push_front(t2_);
  auto const detached1 = fork1->DetachForkWithCopiedBegin();
  EXPECT_TRUE(detached1->is_root());
  auto times = Times(detached1.get());
  EXPECT_THAT(times, ElementsAre(t2_, t4_));
  times = Times(fork2);
  EXPECT_THAT(times, ElementsAre(t1_, t2_));
  times = Times(fork3);
  EXPECT_THAT(times, ElementsAre(t2_));

  fork2->push_front(t2_);
  auto const detached2 = fork2->DetachForkWithCopiedBegin();
  EXPECT_TRUE(detached2->is_root());
  times = Times(detached2.get());
  EXPECT_THAT(times, ElementsAre(t2_));

  fork3->push_front(t2_);
  auto const detached3 = fork3->DetachForkWithCopiedBegin();
  EXPECT_TRUE(detached3->is_root());
  times = Times(detached3.get());
  EXPECT_THAT(times, ElementsAre(t2_));
}

TEST_F(ForkableDeathTest, DeleteAllForksAfterError) {
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    trajectory_.push_back(t2_);
    not_null<FakeTrajectory*> const fork =
        trajectory_.NewFork(trajectory_.timeline_find(t2_));
    fork->DeleteAllForksAfter(t1_);
  }, "before the fork time");
}

TEST_F(ForkableTest, DeleteAllForksAfterSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork =
      trajectory_.NewFork(trajectory_.timeline_find(t2_));
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

TEST_F(ForkableDeathTest, CheckNoForksBeforeError) {
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    not_null<FakeTrajectory*> const fork =
        trajectory_.NewFork(trajectory_.timeline_find(t1_));
    fork->CheckNoForksBefore(t1_);
  }, "nonroot");
  EXPECT_DEATH({
    trajectory_.push_back(t1_);
    trajectory_.push_back(t2_);
    trajectory_.push_back(t3_);
    not_null<FakeTrajectory*> const fork =
        trajectory_.NewFork(trajectory_.timeline_find(t2_));
    fork->push_back(t4_);
    trajectory_.CheckNoForksBefore(t3_);
  }, "found 1 fork");
}

TEST_F(ForkableTest, CheckNoForksBeforeSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork =
      trajectory_.NewFork(trajectory_.timeline_find(t2_));
  fork->push_back(t4_);

  trajectory_.CheckNoForksBefore(t1_ + (t2_ - t1_) / 2);
  auto times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  times = Times(fork);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t4_));

  trajectory_.CheckNoForksBefore(t2_);
  times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
}

TEST_F(ForkableDeathTest, IteratorDecrementError) {
  EXPECT_DEATH({
    auto it = trajectory_.end();
    --it;
  }, "parent_.*non NULL");
}

TEST_F(ForkableTest, IteratorDecrementNoForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  auto it = trajectory_.end();
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
  auto fork = trajectory_.NewFork(trajectory_.timeline_find(t1_));
  trajectory_.push_back(t4_);
  fork->push_back(t3_);
  auto it = fork->end();
  --it;
  EXPECT_EQ(t3_, *it.current());
  --it;
  EXPECT_EQ(t1_, *it.current());
}

TEST_F(ForkableTest, IteratorDecrementMultipleForksSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  auto fork1 = trajectory_.NewFork(trajectory_.timeline_find(t2_));
  auto fork2 = fork1->NewFork(fork1->timeline_find(t2_));
  auto fork3 = fork2->NewFork(fork2->timeline_find(t2_));
  fork2->push_back(t3_);
  auto it = fork3->end();
  --it;
  EXPECT_EQ(t2_, *it.current());
  --it;
  EXPECT_EQ(t1_, *it.current());
  EXPECT_EQ(it, fork3->begin());
}

TEST_F(ForkableDeathTest, IteratorIncrementError) {
  EXPECT_DEATH({
    auto it = trajectory_.begin();
    ++it;
  }, "current.*!=.*end");
}

TEST_F(ForkableTest, IteratorIncrementNoForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  auto it = trajectory_.begin();
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(t3_, *it.current());
}

TEST_F(ForkableTest, IteratorIncrementForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  auto fork = trajectory_.NewFork(trajectory_.timeline_find(t1_));
  trajectory_.push_back(t4_);
  fork->push_back(t3_);
  auto it = fork->begin();
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t3_, *it.current());
  ++it;
  EXPECT_EQ(it, fork->end());
}

TEST_F(ForkableTest, IteratorIncrementMultipleForksSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  auto fork1 = trajectory_.NewFork(trajectory_.timeline_find(t2_));
  auto fork2 = fork1->NewFork(fork1->timeline_find(t2_));
  auto fork3 = fork2->NewFork(fork2->timeline_find(t2_));
  auto it = fork3->begin();
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(it, fork3->end());
  fork3->push_back(t3_);
  --it;
  EXPECT_EQ(t3_, *it.current());
  it = fork3->begin();
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(t3_, *it.current());
  ++it;
  EXPECT_EQ(it, fork3->end());
}

#if !defined(_DEBUG)
TEST_F(ForkableTest, IteratorEndEquality) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  auto fork1 = trajectory_.NewFork(trajectory_.timeline_find(t1_));
  auto fork2 = trajectory_.NewFork(trajectory_.timeline_find(t2_));
  auto it1 = fork1->end();
  auto it2 = fork2->end();
  EXPECT_NE(it1, it2);
}
#endif

TEST_F(ForkableTest, Root) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork =
      trajectory_.NewFork(trajectory_.timeline_find(t2_));
  EXPECT_TRUE(trajectory_.is_root());
  EXPECT_FALSE(fork->is_root());
  EXPECT_EQ(&trajectory_, trajectory_.root());
  EXPECT_EQ(&trajectory_, fork->root());
  EXPECT_EQ(t2_, *fork->Fork().current());
}

TEST_F(ForkableTest, IteratorBeginSuccess) {
  auto it = trajectory_.begin();
  EXPECT_EQ(it, trajectory_.end());

  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);

  it = trajectory_.begin();
  EXPECT_NE(it, trajectory_.end());
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(t3_, *it.current());
  ++it;
  EXPECT_EQ(it, trajectory_.end());

  not_null<FakeTrajectory*> const fork =
      trajectory_.NewFork(trajectory_.timeline_find(t2_));
  fork->push_back(t4_);

  it = fork->begin();
  EXPECT_NE(it, fork->end());
  EXPECT_EQ(t1_, *it.current());
  ++it;
  EXPECT_EQ(t2_, *it.current());
  ++it;
  EXPECT_EQ(t4_, *it.current());
  ++it;
  EXPECT_EQ(it, fork->end());
}

TEST_F(ForkableTest, IteratorFindSuccess) {
  auto it = trajectory_.Find(t0_);
  EXPECT_EQ(it, trajectory_.end());

  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);

  it = trajectory_.Find(t0_);
  EXPECT_EQ(it, trajectory_.end());
  it = trajectory_.Find(t1_);
  EXPECT_NE(it, trajectory_.end());
  EXPECT_EQ(t1_, *it.current());
  it = trajectory_.Find(t2_);
  EXPECT_EQ(t2_, *it.current());
  it = trajectory_.Find(t4_);
  EXPECT_EQ(it, trajectory_.end());

  not_null<FakeTrajectory*> const fork =
      trajectory_.NewFork(trajectory_.timeline_find(t2_));
  fork->push_back(t4_);

  it = fork->Find(t0_);
  EXPECT_EQ(it, fork->end());
  it = fork->Find(t1_);
  EXPECT_NE(it, fork->end());
  EXPECT_EQ(t1_, *it.current());
  it = fork->Find(t2_);
  EXPECT_EQ(t2_, *it.current());
  it = fork->Find(t4_);
  EXPECT_EQ(t4_, *it.current());
  it = fork->Find(t4_ + 1 * Second);
  EXPECT_EQ(it, fork->end());
}

TEST_F(ForkableTest, IteratorLowerBoundSuccess) {
  auto it = trajectory_.LowerBound(t0_);
  EXPECT_EQ(it, trajectory_.end());

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
  EXPECT_EQ(it, trajectory_.end());

  not_null<FakeTrajectory*> const fork1 =
      trajectory_.NewFork(trajectory_.timeline_find(t2_));
  fork1->push_back(t4_);

  it = fork1->LowerBound(t0_);
  EXPECT_EQ(t1_, *it.current());
  it = fork1->LowerBound(t1_);
  EXPECT_NE(it, fork1->end());
  EXPECT_EQ(t1_, *it.current());
  it = fork1->LowerBound(t2_);
  EXPECT_EQ(t2_, *it.current());
  // On the |fork| there is no point at |t3_| so we must land on |t4_|.
  it = fork1->LowerBound(t3_);
  EXPECT_EQ(t4_, *it.current());
  it = fork1->LowerBound(t3_ + 1 * Second);
  EXPECT_EQ(t4_, *it.current());
  it = fork1->LowerBound(t4_);
  EXPECT_EQ(t4_, *it.current());
  it = fork1->LowerBound(t4_ + 1 * Second);
  EXPECT_EQ(it, fork1->end());

  not_null<FakeTrajectory*> const fork2 =
      fork1->NewFork(fork1->timeline_find(t2_));
  fork2->push_back(t5_);
  it = fork2->LowerBound(t4_ - 1 * Second);
  EXPECT_EQ(t5_, *it.current());
}

TEST_F(ForkableTest, IteratorLowerBoundInterestingTimeline) {
  // Test for the "interesting timeline" branch of Iterator::LowerBound.
  trajectory_.push_back(t0_);

  not_null<FakeTrajectory*> const fork =
      trajectory_.NewFork(trajectory_.timeline_find(t0_));
  fork->push_back(t2_);
  not_null<FakeTrajectory*> const fork2 =
      fork->NewFork(fork->timeline_find(t2_));
  fork2->push_back(t3_);

  auto const lower_bound_it = fork2->LowerBound(t1_);
  auto const it1 = ++(fork2->begin());
  EXPECT_EQ(*lower_bound_it, t2_);
  EXPECT_EQ(*it1, t2_);
  EXPECT_EQ(lower_bound_it, it1);
}

TEST_F(ForkableTest, FrontBack) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);

  EXPECT_EQ(t1_, trajectory_.front());
  EXPECT_EQ(t3_, trajectory_.back());

  not_null<FakeTrajectory*> const fork1 =
      trajectory_.NewFork(trajectory_.timeline_find(t2_));
  fork1->push_back(t4_);

  EXPECT_EQ(t1_, fork1->front());
  EXPECT_EQ(t4_, fork1->back());

  not_null<FakeTrajectory*> const fork2 =
      fork1->NewFork(fork1->timeline_find(t2_));
  fork2->push_back(t5_);

  EXPECT_EQ(t1_, fork2->front());
  EXPECT_EQ(t5_, fork2->back());
}

}  // namespace internal_forkable
}  // namespace physics
}  // namespace principia
