#include "trajectory.hpp"

#include "body.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

using principia::geometry::Vector;
using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::Speed;
using principia::quantities::Time;
using principia::quantities::SIUnit;
using principia::si::Metre;
using principia::si::Second;
using testing::ElementsAre;
using testing::Pair;

namespace principia {
namespace physics {

class World;

class TrajectoryTest : public testing::Test {
 protected:
  void SetUp() override {
    q1_ = Vector<Length, World>({1 * Metre, 2 * Metre, 3 * Metre});
    q2_ = Vector<Length, World>({11 * Metre, 12 * Metre, 13 * Metre});
    q3_ = Vector<Length, World>({21 * Metre, 22 * Metre, 23 * Metre});
    q4_ = Vector<Length, World>({31 * Metre, 32 * Metre, 33 * Metre});
    p1_ = Vector<Speed, World>({4 * Metre / Second,
                                5 * Metre / Second,
                                6 * Metre / Second});
    p2_ = Vector<Speed, World>({14 * Metre / Second,
                                15 * Metre / Second,
                                16 * Metre / Second});
    p3_ = Vector<Speed, World>({24 * Metre / Second,
                                25 * Metre / Second,
                                26 * Metre / Second});
    p4_ = Vector<Speed, World>({34 * Metre / Second,
                                35 * Metre / Second,
                                36 * Metre / Second});
    t1_ = 7 * Second;
    t2_ = 17 * Second;
    t3_ = 27 * Second;
    t4_ = 37 * Second;

    body_.reset(new Body(SIUnit<Mass>()));
    trajectory_.reset(new Trajectory<World>(body_.get()));
  }

  Vector<Length, World> q1_, q2_, q3_, q4_;
  Vector<Speed, World> p1_, p2_, p3_, p4_;
  Time t1_, t2_, t3_, t4_;
  std::unique_ptr<Body> body_;
  std::unique_ptr<Trajectory<World>> trajectory_;
};

typedef TrajectoryTest TrajectoryDeathTest;

TEST_F(TrajectoryDeathTest, AppendError) {
  // TODO(phl): Find out how to match the messages here.
  EXPECT_DEATH({
    trajectory_->Append(q2_, p2_, t2_);
    trajectory_->Append(q1_, p1_, t1_);
  }, "");
  EXPECT_DEATH({
    trajectory_->Append(q1_, p1_, t1_);
    trajectory_->Append(q1_, p1_, t1_);
  }, "");
}

TEST_F(TrajectoryTest, AppendSuccess) {
  trajectory_->Append(q1_, p1_, t1_);
  trajectory_->Append(q2_, p2_, t2_);
  trajectory_->Append(q3_, p3_, t3_);
  std::map<Time, Vector<Length, World>> const positions =
      trajectory_->Positions();
  std::map<Time, Vector<Speed, World>> const velocities =
      trajectory_->Velocities();
  std::list<Time> const times = trajectory_->Times();
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(body_.get(), trajectory_->body());
}

TEST_F(TrajectoryTest, ForkError) {
  EXPECT_DEATH({
    trajectory_->Append(q1_, p1_, t1_);
    trajectory_->Append(q3_, p3_, t3_);
    Trajectory<World>* fork = trajectory_->Fork(t2_);
  }, "");
}

TEST_F(TrajectoryTest, ForkSuccess) {
  trajectory_->Append(q1_, p1_, t1_);
  trajectory_->Append(q2_, p2_, t2_);
  trajectory_->Append(q3_, p3_, t3_);
  Trajectory<World>* fork = trajectory_->Fork(t2_);
  fork->Append(q4_, p4_, t4_);
  std::map<Time, Vector<Length, World>> positions = trajectory_->Positions();
  std::map<Time, Vector<Speed, World>> velocities = trajectory_->Velocities();
  std::list<Time> times = trajectory_->Times();
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(body_.get(), fork->body());
  positions = fork->Positions();
  velocities = fork->Velocities();
  times = fork->Times();
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_),
                          Pair(t3_, q3_), Pair(t4_, q4_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_),
                          Pair(t3_, p3_), Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_, t4_));
  EXPECT_EQ(body_.get(), fork->body());
}

TEST_F(TrajectoryTest, Last) {
  trajectory_->Append(q1_, p1_, t1_);
  trajectory_->Append(q2_, p2_, t2_);
  trajectory_->Append(q3_, p3_, t3_);
  EXPECT_EQ(q3_, trajectory_->last_position());
  EXPECT_EQ(p3_, trajectory_->last_velocity());
  EXPECT_EQ(t3_, trajectory_->last_time());
}

TEST_F(TrajectoryTest, Root) {
  trajectory_->Append(q1_, p1_, t1_);
  trajectory_->Append(q2_, p2_, t2_);
  trajectory_->Append(q3_, p3_, t3_);
  Trajectory<World>* fork = trajectory_->Fork(t2_);
  EXPECT_TRUE(trajectory_->is_root());
  EXPECT_FALSE(fork->is_root());
  EXPECT_EQ(trajectory_.get(), trajectory_->root());
  EXPECT_EQ(trajectory_.get(), fork->root());
  EXPECT_EQ(nullptr, trajectory_->fork_time());
  EXPECT_EQ(t2_, *fork->fork_time());
}

TEST_F(TrajectoryTest, ForgetBeforeError) {
  EXPECT_DEATH({
    trajectory_->Append(q1_, p1_, t1_);
    Trajectory<World>* fork = trajectory_->Fork(t1_);
    fork->ForgetBefore(t1_);
  }, "");
  EXPECT_DEATH({
    trajectory_->Append(q1_, p1_, t1_);
    Trajectory<World>* fork = trajectory_->Fork(t1_);
    fork->ForgetBefore(t2_);
  }, "");
}

TEST_F(TrajectoryTest, ForgetBeforeSuccess) {
  trajectory_->Append(q1_, p1_, t1_);
  trajectory_->Append(q2_, p2_, t2_);
  trajectory_->Append(q3_, p3_, t3_);
  Trajectory<World>* fork = trajectory_->Fork(t2_);
  fork->Append(q4_, p4_, t4_);

  trajectory_->ForgetBefore(t1_);
  std::map<Time, Vector<Length, World>> positions = trajectory_->Positions();
  std::map<Time, Vector<Speed, World>> velocities = trajectory_->Velocities();
  std::list<Time> times = trajectory_->Times();
  EXPECT_THAT(positions, ElementsAre(Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t2_, t3_));
  positions = fork->Positions();
  velocities = fork->Velocities();
  times = fork->Times();
  EXPECT_THAT(positions,
              ElementsAre(Pair(t2_, q2_), Pair(t3_, q3_), Pair(t4_, q4_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t2_, p2_), Pair(t3_, p3_), Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t2_, t3_, t4_));

  trajectory_->ForgetBefore(t2_);
  positions = trajectory_->Positions();
  velocities = trajectory_->Velocities();
  times = trajectory_->Times();
  EXPECT_THAT(positions, ElementsAre(Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t3_));
  // Don't use fork, it is dangling.
}

}  // namespace physics
}  // namespace principia
