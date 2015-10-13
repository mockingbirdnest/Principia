#include "physics/discrete_trajectory.hpp"

#include <functional>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::Frame;
using geometry::Instant;
using geometry::Point;
using geometry::R3Element;
using geometry::Vector;
using quantities::Length;
using quantities::Speed;
using quantities::SIUnit;
using quantities::si::Metre;
using quantities::si::Second;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Ref;

// Note that we cannot have a |using ::testing::Pair| here as it would conflict
// with |principia::geometry::Pair|.

namespace physics {

class DiscreteTrajectoryTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  DiscreteTrajectoryTest()
      : q1_(Position<World>(
            Vector<Length, World>({1 * Metre, 2 * Metre, 3 * Metre}))),
        q2_(Position<World>(
            Vector<Length, World>({11 * Metre, 12 * Metre, 13 * Metre}))),
        q3_(Position<World>(
            Vector<Length, World>({21 * Metre, 22 * Metre, 23 * Metre}))),
        q4_(Position<World>(
            Vector<Length, World>({31 * Metre, 32 * Metre, 33 * Metre}))),
        p1_(Velocity<World>({4 * Metre / Second,
                               5 * Metre / Second,
                               6 * Metre / Second})),
        p2_(Velocity<World>({14 * Metre / Second,
                               15 * Metre / Second,
                               16 * Metre / Second})),
        p3_(Velocity<World>({24 * Metre / Second,
                               25 * Metre / Second,
                               26 * Metre / Second})),
        p4_(Velocity<World>({34 * Metre / Second,
                               35 * Metre / Second,
                               36 * Metre / Second})),
        d1_(DegreesOfFreedom<World>(q1_, p1_)),
        d2_(DegreesOfFreedom<World>(q2_, p2_)),
        d3_(DegreesOfFreedom<World>(q3_, p3_)),
        d4_(DegreesOfFreedom<World>(q4_, p4_)) {
    t0_ = Instant(0 * Second);
    t1_ = t0_ + 7 * Second;
    t2_ = t0_ + 17 * Second;
    t3_ = t0_ + 27 * Second;
    t4_ = t0_ + 37 * Second;

    massive_trajectory_ = std::make_unique<DiscreteTrajectory<World>>();
    massless_trajectory_ = std::make_unique<DiscreteTrajectory<World>>();

    transform_ = [](
        Instant const& t,
        DegreesOfFreedom<World> const& from_degrees_of_freedom,
        DiscreteTrajectory<World> const* actual_trajectory,
        DiscreteTrajectory<World> const* expected_trajectory) ->
        DegreesOfFreedom<World> {
      CHECK_EQ(expected_trajectory, actual_trajectory);
      return {Position<World>(
                  2 * (from_degrees_of_freedom.position() -
                       Position<World>(Vector<Length, World>(
                           {43 * Metre, 42 * Metre, 41 * Metre})))),
              3 * from_degrees_of_freedom.velocity()};
    };
    massive_transform_ =
        std::bind(transform_, _1, _2, _3, massive_trajectory_.get());
    massless_transform_ =
        std::bind(transform_, _1, _2, _3, massless_trajectory_.get());
  }

  Position<World> q1_, q2_, q3_, q4_;
  Velocity<World> p1_, p2_, p3_, p4_;
  DegreesOfFreedom<World> d1_, d2_, d3_, d4_;
  Instant t0_, t1_, t2_, t3_, t4_;
  std::unique_ptr<DiscreteTrajectory<World>> massive_trajectory_;
  std::unique_ptr<DiscreteTrajectory<World>> massless_trajectory_;
  std::function<DegreesOfFreedom<World>(
      Instant const&,
      DegreesOfFreedom<World> const&,
      DiscreteTrajectory<World> const*,
      DiscreteTrajectory<World> const*)> transform_;
  DiscreteTrajectory<World>::Transform<World> massive_transform_;
  DiscreteTrajectory<World>::Transform<World> massless_transform_;
};

using DiscreteTrajectoryDeathTest = DiscreteTrajectoryTest;

TEST_F(DiscreteTrajectoryTest, Destruction) {
  int i = 1;
  {
    DiscreteTrajectory<World> massive_trajectory;
  }
  EXPECT_EQ(1, i);
  {
    DiscreteTrajectory<World> massive_trajectory;
    massive_trajectory.set_on_destroy(
        [&i](not_null<DiscreteTrajectory<World>const*> const) { ++i; });
  }
  EXPECT_EQ(2, i);
}

TEST_F(DiscreteTrajectoryTest, NewForkWithCopySuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);
  std::map<Instant, Position<World>> positions =
      massive_trajectory_->Positions();
  std::map<Instant, Velocity<World>> velocities =
      massive_trajectory_->Velocities();
  std::list<Instant> times = massive_trajectory_->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  positions = fork->Positions();
  velocities = fork->Velocities();
  times = fork->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_),
                                     testing::Pair(t4_, q4_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_),
                                      testing::Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_, t4_));
}

TEST_F(DiscreteTrajectoryTest, NewForkWithCopyAtLast) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork1 =
      massive_trajectory_->NewForkWithCopy(t3_);
  not_null<DiscreteTrajectory<World>*> const fork2 =
      fork1->NewForkWithCopy(fork1->last().time());
  not_null<DiscreteTrajectory<World>*> const fork3 =
      fork2->NewForkWithCopy(fork1->last().time());
  EXPECT_EQ(t3_, massive_trajectory_->last().time());
  EXPECT_EQ(t3_, fork1->last().time());

  std::map<Instant, Position<World>> positions = fork2->Positions();
  std::map<Instant, Velocity<World>> velocities = fork2->Velocities();
  std::list<Instant> times = fork2->Times();
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

  after.clear();
  for (auto it = fork2->on_or_after(t3_); !it.at_end(); ++it) {
    after.push_back(it.time());
  }
  EXPECT_THAT(after, ElementsAre(t3_));

  fork1->Append(t4_, d4_);
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

  fork2->Append(t4_, d4_);
  after.clear();
  for (auto it = fork2->on_or_after(t3_); !it.at_end(); ++it) {
    after.push_back(it.time());
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  fork3->Append(t4_, d4_);
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

TEST_F(DiscreteTrajectoryDeathTest, AppendError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t2_, d2_);
    massive_trajectory_->Append(t1_, d1_);
  }, "out of order");
}

TEST_F(DiscreteTrajectoryTest, AppendAtExistingTime) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t1_, d1_);
}

TEST_F(DiscreteTrajectoryTest, AppendSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  std::map<Instant, Position<World>> const positions =
      massive_trajectory_->Positions();
  std::map<Instant, Velocity<World>> const velocities =
      massive_trajectory_->Velocities();
  std::list<Instant> const times = massive_trajectory_->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
}

TEST_F(DiscreteTrajectoryTest, ForgetAfter) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);

  fork->ForgetAfter(t3_ + (t4_ - t3_) / 2);
  std::map<Instant, Position<World>> positions = fork->Positions();
  std::map<Instant, Velocity<World>> velocities = fork->Velocities();
  std::list<Instant> times = fork->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  fork->ForgetAfter(t2_);
  positions = fork->Positions();
  velocities = fork->Velocities();
  times = fork->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_));
  EXPECT_EQ(q2_, fork->last().degrees_of_freedom().position());
  EXPECT_EQ(p2_, fork->last().degrees_of_freedom().velocity());
  EXPECT_EQ(t2_, fork->last().time());

  positions = massive_trajectory_->Positions();
  velocities = massive_trajectory_->Velocities();
  times = massive_trajectory_->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  massive_trajectory_->ForgetAfter(t1_);
  positions = massive_trajectory_->Positions();
  velocities = massive_trajectory_->Velocities();
  times = massive_trajectory_->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_)));
  EXPECT_THAT(times, ElementsAre(t1_));
  // Don't use fork, it is dangling.
}

TEST_F(DiscreteTrajectoryTest, ForgetBefore) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);

  massive_trajectory_->ForgetBefore(t1_ + (t2_ - t1_) / 2);
  std::map<Instant, Position<World>> positions =
      massive_trajectory_->Positions();
  std::map<Instant, Velocity<World>> velocities =
      massive_trajectory_->Velocities();
  std::list<Instant> times = massive_trajectory_->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t2_, t3_));
  positions = fork->Positions();
  velocities = fork->Velocities();
  times = fork->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_),
                                     testing::Pair(t4_, q4_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_),
                                      testing::Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t2_, t3_, t4_));

  massive_trajectory_->ForgetBefore(t2_);
  positions = massive_trajectory_->Positions();
  velocities = massive_trajectory_->Velocities();
  times = massive_trajectory_->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t3_));
  // Don't use fork, it is dangling.
}

TEST_F(DiscreteTrajectoryTest, PointerSerializationSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork1 =
      massive_trajectory_->NewForkWithCopy(t2_);
  not_null<DiscreteTrajectory<World>*> const fork2 =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork2->Append(t4_, d4_);
  not_null<DiscreteTrajectory<World>*> const fork3 =
      massive_trajectory_->NewForkWithCopy(t3_);
  fork3->Append(t4_, d4_);
  serialization::Trajectory root;
  serialization::Trajectory::Pointer root_it;
  serialization::Trajectory::Pointer fork2_it;
  massive_trajectory_->WriteToMessage(&root);
  massive_trajectory_->WritePointerToMessage(&root_it);
  fork2->WritePointerToMessage(&fork2_it);
  EXPECT_EQ(fork2,
            DiscreteTrajectory<World>::ReadPointerFromMessage(
                fork2_it,
                massive_trajectory_.get()));
  EXPECT_EQ(massive_trajectory_.get(),
            DiscreteTrajectory<World>::ReadPointerFromMessage(
                root_it,
                massive_trajectory_.get()));
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const
      massive_trajectory = DiscreteTrajectory<World>::ReadFromMessage(root);
  EXPECT_EQ(massive_trajectory.get(),
            DiscreteTrajectory<World>::ReadPointerFromMessage(
                root_it, massive_trajectory.get()));
}

TEST_F(DiscreteTrajectoryDeathTest, TrajectorySerializationError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    not_null<DiscreteTrajectory<World>*> const fork =
        massive_trajectory_->NewForkWithCopy(t1_);
    serialization::Trajectory message;
    fork->WriteToMessage(&message);
  }, "is_root");
}

TEST_F(DiscreteTrajectoryTest, TrajectorySerializationSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork1 =
      massive_trajectory_->NewForkWithCopy(t2_);
  not_null<DiscreteTrajectory<World>*> const fork2 =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork2->Append(t4_, d4_);
  not_null<DiscreteTrajectory<World>*> const fork3 =
      massive_trajectory_->NewForkWithCopy(t3_);
  fork3->Append(t4_, d4_);
  serialization::Trajectory message;
  serialization::Trajectory reference_message;
  massive_trajectory_->WriteToMessage(&message);
  massive_trajectory_->WriteToMessage(&reference_message);
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const
      deserialized_trajectory =
          DiscreteTrajectory<World>::ReadFromMessage(message);
  message.Clear();
  deserialized_trajectory->WriteToMessage(&message);
  EXPECT_EQ(reference_message.SerializeAsString(), message.SerializeAsString());
  EXPECT_THAT(message.children_size(), Eq(2));
  EXPECT_THAT(message.timeline_size(), Eq(3));
  EXPECT_THAT(Instant::ReadFromMessage(message.timeline(0).instant()), Eq(t1_));
  EXPECT_THAT(Instant::ReadFromMessage(message.timeline(1).instant()), Eq(t2_));
  EXPECT_THAT(Instant::ReadFromMessage(message.timeline(2).instant()), Eq(t3_));
  EXPECT_THAT(
      DegreesOfFreedom<World>::ReadFromMessage(
          message.timeline(0).degrees_of_freedom()), Eq(d1_));
  EXPECT_THAT(
        DegreesOfFreedom<World>::ReadFromMessage(
            message.timeline(1).degrees_of_freedom()), Eq(d2_));
  EXPECT_THAT(
        DegreesOfFreedom<World>::ReadFromMessage(
            message.timeline(2).degrees_of_freedom()), Eq(d3_));
  EXPECT_THAT(message.children(0).trajectories_size(), Eq(2));
  EXPECT_THAT(message.children(0).trajectories(0).children_size(), Eq(0));
  EXPECT_THAT(message.children(0).trajectories(0).timeline_size(), Eq(1));
  EXPECT_THAT(
      Instant::ReadFromMessage(
          message.children(0).trajectories(0).timeline(0).instant()), Eq(t3_));
  EXPECT_THAT(
      DegreesOfFreedom<World>::ReadFromMessage(
          message.children(0).trajectories(0).timeline(0).degrees_of_freedom()),
      Eq(d3_));
  EXPECT_THAT(message.children(0).trajectories(1).children_size(), Eq(0));
  EXPECT_THAT(message.children(0).trajectories(1).timeline_size(), Eq(2));
  EXPECT_THAT(
      Instant::ReadFromMessage(
          message.children(0).trajectories(1).timeline(0).instant()), Eq(t3_));
  EXPECT_THAT(
      DegreesOfFreedom<World>::ReadFromMessage(
          message.children(0).trajectories(1).timeline(0).degrees_of_freedom()),
      Eq(d3_));
  EXPECT_THAT(
      Instant::ReadFromMessage(
          message.children(0).trajectories(1).timeline(1).instant()), Eq(t4_));
  EXPECT_THAT(
      DegreesOfFreedom<World>::ReadFromMessage(
          message.children(0).trajectories(1).timeline(1).degrees_of_freedom()),
      Eq(d4_));
  EXPECT_THAT(message.children(1).trajectories_size(), Eq(1));
  EXPECT_THAT(message.children(1).trajectories(0).children_size(), Eq(0));
  EXPECT_THAT(message.children(1).trajectories(0).timeline_size(), Eq(1));
  EXPECT_THAT(
      Instant::ReadFromMessage(
          message.children(1).trajectories(0).timeline(0).instant()), Eq(t4_));
  EXPECT_THAT(
      DegreesOfFreedom<World>::ReadFromMessage(
          message.children(1).trajectories(0).timeline(0).degrees_of_freedom()),
      Eq(d4_));
}

TEST_F(DiscreteTrajectoryDeathTest, LastError) {
  EXPECT_DEATH({
    massive_trajectory_->last();
  }, "parent_.*non NULL");
}

TEST_F(DiscreteTrajectoryTest, LastSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  EXPECT_EQ(q3_, massive_trajectory_->last().degrees_of_freedom().position());
  EXPECT_EQ(p3_, massive_trajectory_->last().degrees_of_freedom().velocity());
  EXPECT_EQ(t3_, massive_trajectory_->last().time());
}

TEST_F(DiscreteTrajectoryDeathTest, NativeIteratorError) {
  EXPECT_DEATH({
    DiscreteTrajectory<World>::NativeIterator it = massive_trajectory_->last();
  }, "parent_.*non NULL");
}

TEST_F(DiscreteTrajectoryTest, NativeIteratorSuccess) {
  DiscreteTrajectory<World>::NativeIterator it = massive_trajectory_->first();
  EXPECT_TRUE(it.at_end());

  massless_trajectory_->Append(t1_, d1_);
  massless_trajectory_->Append(t2_, d2_);
  massless_trajectory_->Append(t3_, d3_);

  it = massless_trajectory_->first();
  EXPECT_FALSE(it.at_end());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(d1_, it.degrees_of_freedom());
  ++it;
  EXPECT_EQ(t2_, it.time());
  EXPECT_EQ(d2_, it.degrees_of_freedom());
  ++it;
  EXPECT_EQ(t3_, it.time());
  EXPECT_EQ(d3_, it.degrees_of_freedom());
  ++it;
  EXPECT_TRUE(it.at_end());

  not_null<DiscreteTrajectory<World>*> const fork =
      massless_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);

  it = fork->first();
  EXPECT_FALSE(it.at_end());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(d1_, it.degrees_of_freedom());
  ++it;
  EXPECT_EQ(t2_, it.time());
  EXPECT_EQ(d2_, it.degrees_of_freedom());
  ++it;
  EXPECT_EQ(t3_, it.time());
  EXPECT_EQ(d3_, it.degrees_of_freedom());
  ++it;
  EXPECT_EQ(t4_, it.time());
  EXPECT_EQ(d4_, it.degrees_of_freedom());
  ++it;
  EXPECT_TRUE(it.at_end());
}

TEST_F(DiscreteTrajectoryDeathTest, TransformingIteratorError) {
  EXPECT_DEATH({
    DiscreteTrajectory<World>::TransformingIterator<World> it =
        massive_trajectory_->last_with_transform(massive_transform_);
  }, "parent_.*non NULL");
}

TEST_F(DiscreteTrajectoryTest, TransformingIteratorSuccess) {
  DiscreteTrajectory<World>::TransformingIterator<World> it =
      massive_trajectory_->first_with_transform(massive_transform_);
  EXPECT_TRUE(it.at_end());

  massless_trajectory_->Append(t1_, d1_);
  massless_trajectory_->Append(t2_, d2_);
  massless_trajectory_->Append(t3_, d3_);

  it = massless_trajectory_->first_with_transform(massless_transform_);
  EXPECT_FALSE(it.at_end());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-84 * Metre, -80 * Metre, -76 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p1_, it.degrees_of_freedom().velocity());
  ++it;
  EXPECT_EQ(t2_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-64 * Metre, -60 * Metre, -56 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p2_, it.degrees_of_freedom().velocity());
  ++it;
  EXPECT_EQ(t3_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-44 * Metre, -40 * Metre, -36 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p3_, it.degrees_of_freedom().velocity());
  ++it;
  EXPECT_TRUE(it.at_end());

  not_null<DiscreteTrajectory<World>*> const fork =
      massless_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);
  DiscreteTrajectory<World>::Transform<World> const fork_transform =
      std::bind(transform_, _1, _2, _3, fork);

  it = fork->first_with_transform(fork_transform);
  EXPECT_FALSE(it.at_end());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-84 * Metre, -80 * Metre, -76 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p1_, it.degrees_of_freedom().velocity());
  ++it;
  EXPECT_EQ(t2_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-64 * Metre, -60 * Metre, -56 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p2_, it.degrees_of_freedom().velocity());
  ++it;
  EXPECT_EQ(t3_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-44 * Metre, -40 * Metre, -36 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p3_, it.degrees_of_freedom().velocity());
  ++it;
  EXPECT_EQ(t4_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-24 * Metre, -20 * Metre, -16 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p4_, it.degrees_of_freedom().velocity());
  ++it;
  EXPECT_TRUE(it.at_end());
}

TEST_F(DiscreteTrajectoryTest, NativeIteratorOnOrAfterSuccess) {
  DiscreteTrajectory<World>::NativeIterator it =
      massive_trajectory_->on_or_after(t0_);
  EXPECT_TRUE(it.at_end());

  massless_trajectory_->Append(t1_, d1_);
  massless_trajectory_->Append(t2_, d2_);
  massless_trajectory_->Append(t3_, d3_);

  it = massless_trajectory_->on_or_after(t0_);
  EXPECT_FALSE(it.at_end());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(d1_, it.degrees_of_freedom());
  it = massless_trajectory_->on_or_after(t2_);
  EXPECT_EQ(t2_, it.time());
  EXPECT_EQ(d2_, it.degrees_of_freedom());
  it = massless_trajectory_->on_or_after(t4_);
  EXPECT_TRUE(it.at_end());

  not_null<DiscreteTrajectory<World>*> const fork =
      massless_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);

  it = fork->on_or_after(t0_);
  EXPECT_FALSE(it.at_end());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(d1_, it.degrees_of_freedom());
  it = fork->on_or_after(t2_);
  EXPECT_EQ(t2_, it.time());
  EXPECT_EQ(d2_, it.degrees_of_freedom());
  it = fork->on_or_after(t4_);
  EXPECT_EQ(t4_, it.time());
  EXPECT_EQ(d4_, it.degrees_of_freedom());
  it = fork->on_or_after(t4_ + 1 * Second);
  EXPECT_TRUE(it.at_end());
}

TEST_F(DiscreteTrajectoryTest, TransformingIteratorOnOrAfterSuccess) {
  DiscreteTrajectory<World>::TransformingIterator<World> it =
      massive_trajectory_->on_or_after_with_transform(t0_, massive_transform_);
  EXPECT_TRUE(it.at_end());

  massless_trajectory_->Append(t1_, d1_);
  massless_trajectory_->Append(t2_, d2_);
  massless_trajectory_->Append(t3_, d3_);

  it = massless_trajectory_->on_or_after_with_transform(t0_,
                                                        massless_transform_);
  EXPECT_FALSE(it.at_end());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-84 * Metre, -80 * Metre, -76 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p1_, it.degrees_of_freedom().velocity());
  it = massless_trajectory_->on_or_after_with_transform(t2_,
                                                        massless_transform_);
  EXPECT_EQ(t2_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-64 * Metre, -60 * Metre, -56 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p2_, it.degrees_of_freedom().velocity());
  it = massless_trajectory_->on_or_after_with_transform(t4_,
                                                        massless_transform_);
  EXPECT_TRUE(it.at_end());

  not_null<DiscreteTrajectory<World>*> const fork =
      massless_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);
  DiscreteTrajectory<World>::Transform<World> const fork_transform =
      std::bind(transform_, _1, _2, _3, fork);

  it = fork->on_or_after_with_transform(t0_, fork_transform);
  EXPECT_FALSE(it.at_end());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(t1_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-84 * Metre, -80 * Metre, -76 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p1_, it.degrees_of_freedom().velocity());
  it = fork->on_or_after_with_transform(t2_, fork_transform);
  EXPECT_EQ(t2_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-64 * Metre, -60 * Metre, -56 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p2_, it.degrees_of_freedom().velocity());
  it = fork->on_or_after_with_transform(t4_, fork_transform);
  EXPECT_EQ(t4_, it.time());
  EXPECT_EQ(Position<World>(Vector<Length, World>(
                {-24 * Metre, -20 * Metre, -16 * Metre})),
            it.degrees_of_freedom().position());
  EXPECT_EQ(3 * p4_, it.degrees_of_freedom().velocity());
  it = fork->on_or_after_with_transform(t4_ + 1 * Second, fork_transform);
  EXPECT_TRUE(it.at_end());
}

}  // namespace physics
}  // namespace principia
