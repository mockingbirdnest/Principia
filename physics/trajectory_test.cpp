#include "trajectory.hpp"

#include <list>
#include <map>

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
using testing::Ref;

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
    d1_ = std::make_unique<DegreesOfFreedom<World>>(q1_, p1_);
    d2_ = std::make_unique<DegreesOfFreedom<World>>(q2_, p2_);
    d3_ = std::make_unique<DegreesOfFreedom<World>>(q3_, p3_);
    d4_ = std::make_unique<DegreesOfFreedom<World>>(q4_, p4_);
    t1_ = 7 * Second;
    t2_ = 17 * Second;
    t3_ = 27 * Second;
    t4_ = 37 * Second;

    body_.reset(new Body(SIUnit<Mass>()));
    trajectory_.reset(new Trajectory<World>(*body_));
  }

  Vector<Length, World> q1_, q2_, q3_, q4_;
  Vector<Speed, World> p1_, p2_, p3_, p4_;
  std::unique_ptr<DegreesOfFreedom<World>> d1_, d2_, d3_, d4_;
  Time t1_, t2_, t3_, t4_;
  std::unique_ptr<Body> body_;
  std::unique_ptr<Trajectory<World>> trajectory_;
};

typedef TrajectoryTest TrajectoryDeathTest;

TEST_F(TrajectoryDeathTest, AppendError) {
  EXPECT_DEATH({
    trajectory_->Append(t2_, *d2_);
    trajectory_->Append(t1_, *d1_);
  }, "out of order");
  EXPECT_DEATH({
    trajectory_->Append(t1_, *d1_);
    trajectory_->Append(t1_, *d1_);
  }, "existing time");
}

TEST_F(TrajectoryTest, AppendSuccess) {
  trajectory_->Append(t1_, *d1_);
  trajectory_->Append(t2_, *d2_);
  trajectory_->Append(t3_, *d3_);
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
  EXPECT_THAT(trajectory_->body(), Ref(*body_));
}

TEST_F(TrajectoryDeathTest, ForkError) {
  EXPECT_DEATH({
    trajectory_->Append(t1_, *d1_);
    trajectory_->Append(t3_, *d3_);
    Trajectory<World>* fork = trajectory_->Fork(t2_);
  }, "nonexistent time");
}

TEST_F(TrajectoryTest, ForkSuccess) {
  trajectory_->Append(t1_, *d1_);
  trajectory_->Append(t2_, *d2_);
  trajectory_->Append(t3_, *d3_);
  Trajectory<World>* fork = trajectory_->Fork(t2_);
  fork->Append(t4_, *d4_);
  std::map<Time, Vector<Length, World>> positions = trajectory_->Positions();
  std::map<Time, Vector<Speed, World>> velocities = trajectory_->Velocities();
  std::list<Time> times = trajectory_->Times();
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_THAT(fork->body(), Ref(*body_));
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
  EXPECT_THAT(fork->body(), Ref(*body_));
}

TEST_F(TrajectoryTest, Last) {
  trajectory_->Append(t1_, *d1_);
  trajectory_->Append(t2_, *d2_);
  trajectory_->Append(t3_, *d3_);
  EXPECT_EQ(q3_, trajectory_->last_position());
  EXPECT_EQ(p3_, trajectory_->last_velocity());
  EXPECT_EQ(t3_, trajectory_->last_time());
}

TEST_F(TrajectoryTest, Root) {
  trajectory_->Append(t1_, *d1_);
  trajectory_->Append(t2_, *d2_);
  trajectory_->Append(t3_, *d3_);
  Trajectory<World>* fork = trajectory_->Fork(t2_);
  EXPECT_TRUE(trajectory_->is_root());
  EXPECT_FALSE(fork->is_root());
  EXPECT_EQ(trajectory_.get(), trajectory_->root());
  EXPECT_EQ(trajectory_.get(), fork->root());
  EXPECT_EQ(nullptr, trajectory_->fork_time());
  EXPECT_EQ(t2_, *fork->fork_time());
}

TEST_F(TrajectoryDeathTest, ForgetAfterError) {
  EXPECT_DEATH({
    trajectory_->Append(t1_, *d1_);
    trajectory_->ForgetAfter(t2_);
  }, "nonexistent time.* root");
  EXPECT_DEATH({
    trajectory_->Append(t1_, *d1_);
    Trajectory<World>* fork = trajectory_->Fork(t1_);
    fork->ForgetAfter(t2_);
  }, "nonexistent time.* nonroot");
}

TEST_F(TrajectoryTest, ForgetAfterSuccess) {
  trajectory_->Append(t1_, *d1_);
  trajectory_->Append(t2_, *d2_);
  trajectory_->Append(t3_, *d3_);
  Trajectory<World>* fork = trajectory_->Fork(t2_);
  fork->Append(t4_, *d4_);

  fork->ForgetAfter(t3_);
  std::map<Time, Vector<Length, World>> positions = fork->Positions();
  std::map<Time, Vector<Speed, World>> velocities = fork->Velocities();
  std::list<Time> times = fork->Times();
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  fork->ForgetAfter(t2_);
  positions = fork->Positions();
  velocities = fork->Velocities();
  times = fork->Times();
  EXPECT_THAT(positions, ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_));

  positions = trajectory_->Positions();
  velocities = trajectory_->Velocities();
  times = trajectory_->Times();
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  trajectory_->ForgetAfter(t1_);
  positions = trajectory_->Positions();
  velocities = trajectory_->Velocities();
  times = trajectory_->Times();
  EXPECT_THAT(positions, ElementsAre(Pair(t1_, q1_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t1_, p1_)));
  EXPECT_THAT(times, ElementsAre(t1_));
  // Don't use fork, it is dangling.
}

TEST_F(TrajectoryDeathTest, ForgetBeforeError) {
  EXPECT_DEATH({
    trajectory_->Append(t1_, *d1_);
    Trajectory<World>* fork = trajectory_->Fork(t1_);
    fork->ForgetBefore(t1_);
  }, "nonroot");
  EXPECT_DEATH({
    trajectory_->Append(t1_, *d1_);
    trajectory_->ForgetBefore(t2_);
  }, "nonexistent time");
}

TEST_F(TrajectoryTest, ForgetBeforeSuccess) {
  trajectory_->Append(t1_, *d1_);
  trajectory_->Append(t2_, *d2_);
  trajectory_->Append(t3_, *d3_);
  Trajectory<World>* fork = trajectory_->Fork(t2_);
  fork->Append(t4_, *d4_);

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
