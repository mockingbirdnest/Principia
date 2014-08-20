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

namespace principia {
namespace physics {

class World;

class TrajectoryTest : public testing::Test {
 protected:
  void SetUp() override {
    q1_ = Vector<Length, World>({1 * Metre, 2 * Metre, 3 * Metre});
    q2_ = Vector<Length, World>({11 * Metre, 12 * Metre, 13 * Metre});
    q3_ = Vector<Length, World>({21 * Metre, 22 * Metre, 23 * Metre});
    p1_ = Vector<Speed, World>({4 * Metre / Second,
                                5 * Metre / Second,
                                6 * Metre / Second});
    p2_ = Vector<Speed, World>({14 * Metre / Second,
                                15 * Metre / Second,
                                16 * Metre / Second});
    p3_ = Vector<Speed, World>({24 * Metre / Second,
                                25 * Metre / Second,
                                26 * Metre / Second});
    t1_ = 7 * Second;
    t2_ = 17 * Second;
    t3_ = 27 * Second;

    body_.reset(new Body(SIUnit<Mass>()));
    trajectory_.reset(new Trajectory<World>(body_.get()));
  }

  Vector<Length, World> q1_, q2_, q3_;
  Vector<Speed, World> p1_, p2_, p3_;
  Time t1_, t2_, t3_;
  std::unique_ptr<Body> body_;
  std::unique_ptr<Trajectory<World>> trajectory_;
};

TEST_F(TrajectoryTest, AppendError) {
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
  std::vector<Vector<Length, World>> const positions = trajectory_->positions();
  std::vector<Vector<Speed, World>> const velocities =
      trajectory_->velocities();
  std::vector<Time> const times = trajectory_->times();
  EXPECT_THAT(positions, ElementsAre(q1_, q2_, q3_));
  EXPECT_THAT(velocities, ElementsAre(p1_, p2_, p3_));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(body_.get(), trajectory_->body());
}

TEST_F(TrajectoryTest, Fork) {
  trajectory_->Append(q1_, p1_, t1_);
  trajectory_->Append(q2_, p2_, t2_);
  trajectory_->Append(q3_, p3_, t3_);
  Trajectory<World>* fork = trajectory_->Fork(t2_);
  std::vector<Vector<Length, World>> const positions = fork->positions();
  std::vector<Vector<Speed, World>> const velocities = fork->velocities();
  std::vector<Time> const times = fork->times();
  EXPECT_THAT(positions, ElementsAre(q1_, q2_, q3_));
  EXPECT_THAT(velocities, ElementsAre(p1_, p2_, p3_));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(body_.get(), fork->body());
}

}  // namespace physics
}  // namespace principia
