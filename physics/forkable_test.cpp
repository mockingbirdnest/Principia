#include "physics/forkable.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {

using geometry::Instant;
using si::Second;
using ::testing::ElementsAre;

namespace physics {

class FakeTrajectory : public Forkable<FakeTrajectory,
                                       std::vector<Instant>::const_iterator> {
 public:
  FakeTrajectory();

  std::vector<Instant>::const_iterator begin() const;
  std::vector<Instant>::const_iterator end() const;

  void push_back(Instant const& time);

 protected:
  not_null<FakeTrajectory*> that() override;
  not_null<FakeTrajectory const*> that() const override;

  TimelineConstIterator timeline_begin() const override;
  TimelineConstIterator timeline_end() const override;
  TimelineConstIterator timeline_find(Instant const& time) const override;
  void timeline_insert(TimelineConstIterator begin,
                       TimelineConstIterator end) override;
  bool timeline_empty() const override;

 private:
  std::vector<Instant> timeline_;

  template<typename Tr4jectory, typename TimelineConstIterator_>
  friend class Forkable;
  template<typename Tr4jectory, typename TimelineConstIterator_>
  friend class Forkable<Tr4jectory, TimelineConstIterator_>::Iterator;
};

FakeTrajectory::FakeTrajectory()
    : Forkable<FakeTrajectory,
               std::vector<Instant>::const_iterator>() {}

std::vector<Instant>::const_iterator FakeTrajectory::begin() const {
  return timeline_.begin();
}

std::vector<Instant>::const_iterator FakeTrajectory::end() const {
  return timeline_.end();
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

  static std::vector<Instant> Times(
      not_null<FakeTrajectory const*> const trajectory) {
    //TODO(phl): Not very nice.
    std::vector<Instant> times;
    FakeTrajectory::Iterator it =
        FakeTrajectory::Iterator::New(trajectory,
                                      trajectory->root(),
                                      trajectory->root()->begin());
    for (; it.current() != trajectory->end(); ++it) {
      times.push_back(*it.current());
    }
    return times;
  }

  static Instant const& LastTime(
      not_null<FakeTrajectory const*> const trajectory) {
    FakeTrajectory::Iterator it =
        FakeTrajectory::Iterator::New(trajectory,
                                      trajectory,
                                      --trajectory->end());
    return *it.current();
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

TEST_F(ForkableDeathTest, ForkSuccess) {
  trajectory_.push_back(t1_);
  trajectory_.push_back(t2_);
  trajectory_.push_back(t3_);
  not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
  fork->push_back(t4_);
  auto times = Times(&trajectory_);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  times = Times(fork);
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_, t4_));
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
  EXPECT_EQ(t3_, *fork2->ForkTime());

  std::vector<Instant> after;
  for (auto it = fork3->on_or_after(t3_); !it.at_end(); ++it) {
    after.push_back(it.time());
  }
  EXPECT_THAT(after, ElementsAre(t3_));

  fork2->ForgetAfter(t3_);
  positions = fork2->Positions();
  velocities = fork2->Velocities();
  times = fork2->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork2->last().degrees_of_freedom().position());
  EXPECT_EQ(p3_, fork2->last().degrees_of_freedom().velocity());
  EXPECT_EQ(t3_, fork2->last().time());
  EXPECT_EQ(t3_, *fork2->fork_time());

  after.clear();
  for (auto it = fork2->on_or_after(t3_); !it.at_end(); ++it) {
    after.push_back(it.time());
  }
  EXPECT_THAT(after, ElementsAre(t3_));

  fork1.push_back(t4_);
  positions = fork2->Positions();
  velocities = fork2->Velocities();
  times = fork2->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork2->last().degrees_of_freedom().position());
  EXPECT_EQ(p3_, fork2->last().degrees_of_freedom().velocity());
  EXPECT_EQ(t3_, fork2->last().time());
  EXPECT_EQ(t3_, *fork2->fork_time());

  after.clear();
  for (auto it = fork1->on_or_after(t3_); !it.at_end(); ++it) {
    after.push_back(it.time());
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  positions = fork3->Positions();
  velocities = fork3->Velocities();
  times = fork3->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork3->last().degrees_of_freedom().position());
  EXPECT_EQ(p3_, fork3->last().degrees_of_freedom().velocity());
  EXPECT_EQ(t3_, fork3->last().time());
  EXPECT_EQ(t3_, *fork3->fork_time());

  fork2.push_back(t4_);
  after.clear();
  for (auto it = fork2->on_or_after(t3_); !it.at_end(); ++it) {
    after.push_back(it.time());
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  fork3.push_back(t4_);
  after.clear();
  for (auto it = fork3->on_or_after(t3_); !it.at_end(); ++it) {
    after.push_back(it.time());
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  after.clear();
  for (auto it = fork3->on_or_after(t2_); !it.at_end(); ++it) {
    after.push_back(it.time());
  }
  EXPECT_THAT(after, ElementsAre(t2_, t3_, t4_));
}

//TEST_F(TrajectoryDeathTest, DeleteForkError) {
//  EXPECT_DEATH({
//    trajectory_.push_back(t1_);
//    FakeTrajectory* root = trajectory_.get();
//    trajectory_.DeleteFork(&root);
//  }, "'fork_time'.* non NULL");
//  EXPECT_DEATH({
//    trajectory_.push_back(t1_);
//    FakeTrajectory* fork1 = trajectory_.NewFork(t1_);
//    fork1.push_back(t2_);
//    FakeTrajectory* fork2 = fork1->NewFork(t2_);
//    trajectory_.DeleteFork(&fork2);
//  }, "not a child");
//}
//
//TEST_F(ForkableDeathTest, DeleteForkSuccess) {
//  trajectory_.push_back(t1_);
//  trajectory_.push_back(t2_);
//  trajectory_.push_back(t3_);
//  not_null<FakeTrajectory*> const fork1 = trajectory_.NewFork(t2_);
//  FakeTrajectory* fork2 = trajectory_.NewFork(t2_);
//  fork1.push_back(t4_);
//  trajectory_.DeleteFork(&fork2);
//  EXPECT_EQ(nullptr, fork2);
//  std::map<Instant, Position<World>> positions =
//      trajectory_.Positions();
//  std::map<Instant, Velocity<World>> velocities =
//      trajectory_.Velocities();
//  std::list<Instant> times = trajectory_.Times();
//  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
//                                     testing::Pair(t2_, q2_),
//                                     testing::Pair(t3_, q3_)));
//  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
//                                      testing::Pair(t2_, p2_),
//                                      testing::Pair(t3_, p3_)));
//  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
//  EXPECT_THAT(fork1->body<MassiveBody>(), Eq(&massive_body_));
//  positions = fork1->Positions();
//  velocities = fork1->Velocities();
//  times = fork1->Times();
//  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
//                                     testing::Pair(t2_, q2_),
//                                     testing::Pair(t3_, q3_),
//                                     testing::Pair(t4_, q4_)));
//  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
//                                      testing::Pair(t2_, p2_),
//                                      testing::Pair(t3_, p3_),
//                                      testing::Pair(t4_, p4_)));
//  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_, t4_));
//  EXPECT_THAT(fork1->body<MassiveBody>(), Eq(&massive_body_));
//  trajectory_.reset();
//}
//
//TEST_F(TrajectoryDeathTest, LastError) {
//  EXPECT_DEATH({
//    trajectory_.last();
//  }, "Empty trajectory");
//}
//
//TEST_F(ForkableDeathTest, LastSuccess) {
//  trajectory_.push_back(t1_);
//  trajectory_.push_back(t2_);
//  trajectory_.push_back(t3_);
//  EXPECT_EQ(q3_, trajectory_.last().degrees_of_freedom().position());
//  EXPECT_EQ(p3_, trajectory_.last().degrees_of_freedom().velocity());
//  EXPECT_EQ(t3_, trajectory_.last().time());
//}
//
//TEST_F(ForkableDeathTest, Root) {
//  trajectory_.push_back(t1_);
//  trajectory_.push_back(t2_);
//  trajectory_.push_back(t3_);
//  not_null<FakeTrajectory*> const fork = trajectory_.NewFork(t2_);
//  EXPECT_TRUE(trajectory_.is_root());
//  EXPECT_FALSE(fork->is_root());
//  EXPECT_EQ(trajectory_.get(), trajectory_.root());
//  EXPECT_EQ(trajectory_.get(), fork->root());
//  EXPECT_EQ(nullptr, trajectory_.fork_time());
//  EXPECT_EQ(t2_, *fork->fork_time());
//}
//
//TEST_F(TrajectoryDeathTest, NativeIteratorError) {
//  EXPECT_DEATH({
//    FakeTrajectory::NativeIterator it = trajectory_.last();
//  }, "Empty trajectory");
//  EXPECT_DEATH({
//    FakeTrajectory::NativeIterator it = trajectory_.first();
//    ++it;
//  }, "beyond end");
//}
//
//TEST_F(ForkableDeathTest, NativeIteratorSuccess) {
//  FakeTrajectory::NativeIterator it = trajectory_.first();
//  EXPECT_TRUE(it.at_end());
//
//  massless_trajectory_.push_back(t1_);
//  massless_trajectory_.push_back(t2_);
//  massless_trajectory_.push_back(t3_);
//
//  it = massless_trajectory_.first();
//  EXPECT_FALSE(it.at_end());
//  EXPECT_EQ(t1_, it.time());
//  EXPECT_EQ(d1_, it.degrees_of_freedom());
//  ++it;
//  EXPECT_EQ(t2_, it.time());
//  EXPECT_EQ(d2_, it.degrees_of_freedom());
//  ++it;
//  EXPECT_EQ(t3_, it.time());
//  EXPECT_EQ(d3_, it.degrees_of_freedom());
//  ++it;
//  EXPECT_TRUE(it.at_end());
//
//  not_null<FakeTrajectory*> const fork = massless_trajectory_.NewFork(t2_);
//  fork.push_back(t4_);
//
//  it = fork->first();
//  EXPECT_FALSE(it.at_end());
//  EXPECT_EQ(t1_, it.time());
//  EXPECT_EQ(d1_, it.degrees_of_freedom());
//  ++it;
//  EXPECT_EQ(t2_, it.time());
//  EXPECT_EQ(d2_, it.degrees_of_freedom());
//  ++it;
//  EXPECT_EQ(t3_, it.time());
//  EXPECT_EQ(d3_, it.degrees_of_freedom());
//  ++it;
//  EXPECT_EQ(t4_, it.time());
//  EXPECT_EQ(d4_, it.degrees_of_freedom());
//  ++it;
//  EXPECT_TRUE(it.at_end());
//}
//
//TEST_F(ForkableDeathTest, NativeIteratorOnOrAfterSuccess) {
//  FakeTrajectory::NativeIterator it = trajectory_.on_or_after(t0_);
//  EXPECT_TRUE(it.at_end());
//
//  massless_trajectory_.push_back(t1_);
//  massless_trajectory_.push_back(t2_);
//  massless_trajectory_.push_back(t3_);
//
//  it = massless_trajectory_.on_or_after(t0_);
//  EXPECT_FALSE(it.at_end());
//  EXPECT_EQ(t1_, it.time());
//  EXPECT_EQ(d1_, it.degrees_of_freedom());
//  it = massless_trajectory_.on_or_after(t2_);
//  EXPECT_EQ(t2_, it.time());
//  EXPECT_EQ(d2_, it.degrees_of_freedom());
//  it = massless_trajectory_.on_or_after(t4_);
//  EXPECT_TRUE(it.at_end());
//
//  not_null<FakeTrajectory*> const fork = massless_trajectory_.NewFork(t2_);
//  fork.push_back(t4_);
//
//  it = fork->on_or_after(t0_);
//  EXPECT_FALSE(it.at_end());
//  EXPECT_EQ(t1_, it.time());
//  EXPECT_EQ(d1_, it.degrees_of_freedom());
//  it = fork->on_or_after(t2_);
//  EXPECT_EQ(t2_, it.time());
//  EXPECT_EQ(d2_, it.degrees_of_freedom());
//  it = fork->on_or_after(t4_);
//  EXPECT_EQ(t4_, it.time());
//  EXPECT_EQ(d4_, it.degrees_of_freedom());
//  it = fork->on_or_after(t4_ + 1 * Second);
//  EXPECT_TRUE(it.at_end());
//}

}  // namespace physics
}  // namespace principia
