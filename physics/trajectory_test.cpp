#include "trajectory.hpp"

#include <functional>
#include <list>
#include <map>
#include <string>

#include "body.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/frame.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

using principia::geometry::Instant;
using principia::geometry::Point;
using principia::geometry::R3Element;
using principia::geometry::Vector;
using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::Speed;
using principia::quantities::SIUnit;
using principia::si::Metre;
using principia::si::Second;
using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using testing::ElementsAre;
using testing::Eq;
using testing::Ref;

// Note that we cannot have a |using testing::Pair| here as it would conflict
// with |principia::geometry::Pair|.

namespace principia {
namespace physics {

class TrajectoryTest : public testing::Test {
 protected:
  enum class Tag {
    kWorld,
    kOtherWorld,
  };

  using World = Frame<Tag, Tag::kWorld, true>;

  TrajectoryTest()
      : massive_body_(MassiveBody(1 * SIUnit<Mass>())),
        q1_(Position<World>(
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

    massive_trajectory_ = std::make_unique<Trajectory<World>>(&massive_body_);
    massless_trajectory_ = std::make_unique<Trajectory<World>>(&massless_body_);

    transform_ = [](
        Instant const& t,
        DegreesOfFreedom<World> const& from_degrees_of_freedom,
        Trajectory<World> const* actual_trajectory,
        Trajectory<World> const* expected_trajectory) ->
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

  MassiveBody massive_body_;
  MasslessBody massless_body_;
  Position<World> q1_, q2_, q3_, q4_;
  Velocity<World> p1_, p2_, p3_, p4_;
  DegreesOfFreedom<World> d1_, d2_, d3_, d4_;
  Instant t0_, t1_, t2_, t3_, t4_;
  std::unique_ptr<Trajectory<World>> massive_trajectory_;
  std::unique_ptr<Trajectory<World>> massless_trajectory_;
  std::function<DegreesOfFreedom<World>(Instant const&,
                                        DegreesOfFreedom<World> const&,
                                        Trajectory<World> const*,
                                        Trajectory<World> const*)> transform_;
  Trajectory<World>::Transform<World> massive_transform_;
  Trajectory<World>::Transform<World> massless_transform_;
};

using TrajectoryDeathTest = TrajectoryTest;

TEST_F(TrajectoryDeathTest, Construction) {
  using OtherWorld = Frame<Tag, Tag::kOtherWorld, true>;
  EXPECT_DEATH({
    OblateBody<OtherWorld> body(1 * SIUnit<GravitationalParameter>(),
                                1.0 /*j2*/,
                                1 * SIUnit<Length>(),
                                Vector<double, OtherWorld>({0, 1, 0}));
    Trajectory<World> trajectory(&body);
  }, "not in the same frame");
}

TEST_F(TrajectoryDeathTest, AppendError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t2_, d2_);
    massive_trajectory_->Append(t1_, d1_);
  }, "out of order");
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    massive_trajectory_->Append(t1_, d1_);
  }, "existing time");
}

TEST_F(TrajectoryTest, AppendSuccess) {
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
  EXPECT_THAT(massive_trajectory_->body<MassiveBody>(), Eq(&massive_body_));
}

TEST_F(TrajectoryDeathTest, ForkError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    massive_trajectory_->Append(t3_, d3_);
    massive_trajectory_->Fork(t2_);
  }, "nonexistent time");
}

TEST_F(TrajectoryTest, ForkSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<Trajectory<World>*> const fork = massive_trajectory_->Fork(t2_);
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
  EXPECT_THAT(fork->body<MassiveBody>(), Eq(&massive_body_));
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
  EXPECT_THAT(fork->body<MassiveBody>(), Eq(&massive_body_));
}

TEST_F(TrajectoryTest, Serialization) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<Trajectory<World>*> const fork1 = massive_trajectory_->Fork(t1_);
  not_null<Trajectory<World>*> const fork2 = massive_trajectory_->Fork(t1_);
  not_null<Trajectory<World>*> const fork3 = massive_trajectory_->Fork(t3_);
  fork3->Append(t4_, d4_);
  serialization::Trajectory message;
  massive_trajectory_->WriteToMessage(&message);
  EXPECT_THAT(message.children_size(), Eq(2));
  EXPECT_THAT(message.timeline_size(), Eq(3));
  EXPECT_THAT(message.children(0).trajectories_size(), Eq(2));
  EXPECT_THAT(message.children(0).trajectories(0).children_size(), Eq(0));
  EXPECT_THAT(message.children(0).trajectories(1).children_size(), Eq(0));
  EXPECT_THAT(message.children(0).trajectories(0).timeline_size(), Eq(2));
  EXPECT_THAT(message.children(0).trajectories(1).timeline_size(), Eq(2));
  EXPECT_THAT(message.children(1).trajectories_size(), Eq(1));
  EXPECT_THAT(message.children(1).trajectories(0).children_size(), Eq(0));
  EXPECT_THAT(message.children(1).trajectories(0).timeline_size(), Eq(1));
}

TEST_F(TrajectoryDeathTest, DeleteForkError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    Trajectory<World>* root = massive_trajectory_.get();
    massive_trajectory_->DeleteFork(&root);
  }, "'fork_time'.* non NULL");
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    Trajectory<World>* fork1 = massive_trajectory_->Fork(t1_);
    fork1->Append(t2_, d2_);
    Trajectory<World>* fork2 = fork1->Fork(t2_);
    massive_trajectory_->DeleteFork(&fork2);
  }, "not a child");
}

TEST_F(TrajectoryTest, DeleteForkSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<Trajectory<World>*> const fork1 = massive_trajectory_->Fork(t2_);
  Trajectory<World>* fork2 = massive_trajectory_->Fork(t2_);
  fork1->Append(t4_, d4_);
  massive_trajectory_->DeleteFork(&fork2);
  EXPECT_EQ(nullptr, fork2);
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
  EXPECT_THAT(fork1->body<MassiveBody>(), Eq(&massive_body_));
  positions = fork1->Positions();
  velocities = fork1->Velocities();
  times = fork1->Times();
  EXPECT_THAT(positions, ElementsAre(testing::Pair(t1_, q1_),
                                     testing::Pair(t2_, q2_),
                                     testing::Pair(t3_, q3_),
                                     testing::Pair(t4_, q4_)));
  EXPECT_THAT(velocities, ElementsAre(testing::Pair(t1_, p1_),
                                      testing::Pair(t2_, p2_),
                                      testing::Pair(t3_, p3_),
                                      testing::Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_, t4_));
  EXPECT_THAT(fork1->body<MassiveBody>(), Eq(&massive_body_));
  massive_trajectory_.reset();
}

TEST_F(TrajectoryDeathTest, LastError) {
  EXPECT_DEATH({
    massive_trajectory_->last();
  }, "Empty trajectory");
}

TEST_F(TrajectoryTest, LastSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  EXPECT_EQ(q3_, massive_trajectory_->last().degrees_of_freedom().position());
  EXPECT_EQ(p3_, massive_trajectory_->last().degrees_of_freedom().velocity());
  EXPECT_EQ(t3_, massive_trajectory_->last().time());
}

TEST_F(TrajectoryTest, Root) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<Trajectory<World>*> const fork = massive_trajectory_->Fork(t2_);
  EXPECT_TRUE(massive_trajectory_->is_root());
  EXPECT_FALSE(fork->is_root());
  EXPECT_EQ(massive_trajectory_.get(), massive_trajectory_->root());
  EXPECT_EQ(massive_trajectory_.get(), fork->root());
  EXPECT_EQ(nullptr, massive_trajectory_->fork_time());
  EXPECT_EQ(t2_, *fork->fork_time());
}

TEST_F(TrajectoryDeathTest, ForgetAfterError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    massive_trajectory_->ForgetAfter(t2_);
  }, "nonexistent time.* root");
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    not_null<Trajectory<World>*> const fork = massive_trajectory_->Fork(t1_);
    fork->ForgetAfter(t2_);
  }, "nonexistent time.* nonroot");
}

TEST_F(TrajectoryTest, ForgetAfterSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<Trajectory<World>*> const fork = massive_trajectory_->Fork(t2_);
  fork->Append(t4_, d4_);

  fork->ForgetAfter(t3_);
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

TEST_F(TrajectoryDeathTest, ForgetBeforeError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    not_null<Trajectory<World>*> const fork = massive_trajectory_->Fork(t1_);
    fork->ForgetBefore(t1_);
  }, "nonroot");
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    massive_trajectory_->ForgetBefore(t2_);
  }, "nonexistent time");
}

TEST_F(TrajectoryTest, ForgetBeforeSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<Trajectory<World>*> const fork = massive_trajectory_->Fork(t2_);
  fork->Append(t4_, d4_);

  massive_trajectory_->ForgetBefore(t1_);
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

TEST_F(TrajectoryDeathTest, IntrinsicAccelerationError) {
  EXPECT_DEATH({
    massive_trajectory_->set_intrinsic_acceleration(
        [](Instant const& t) { return Vector<Acceleration, World>(); } );
  }, "massive body");
  EXPECT_DEATH({
    massless_trajectory_->set_intrinsic_acceleration(
        [](Instant const& t) { return Vector<Acceleration, World>(); } );
    massless_trajectory_->set_intrinsic_acceleration(
        [](Instant const& t) { return Vector<Acceleration, World>(); } );
  }, "already has.* acceleration");
}

TEST_F(TrajectoryDeathTest, IntrinsicAccelerationSuccess) {
  massless_trajectory_->Append(t1_, d1_);
  massless_trajectory_->Append(t2_, d2_);
  massless_trajectory_->Append(t3_, d3_);
  not_null<Trajectory<World>*> const fork = massless_trajectory_->Fork(t2_);
  fork->Append(t4_, d4_);

  EXPECT_FALSE(massless_trajectory_->has_intrinsic_acceleration());
  massless_trajectory_->set_intrinsic_acceleration(
      [this](Instant const& t) {
        return Vector<Acceleration, World>(
            {1 * SIUnit<Length>() / ((t - t0_) * (t - t0_)),
             2 * SIUnit<Length>() / ((t - t0_) * (t - t0_)),
             3 * SIUnit<Length>() / ((t - t0_) * (t - t0_))});});
  EXPECT_TRUE(massless_trajectory_->has_intrinsic_acceleration());
  EXPECT_THAT(massless_trajectory_->evaluate_intrinsic_acceleration(t1_),
              Eq(Vector<Acceleration, World>(
                  {0.020408163265306122449 * SIUnit<Acceleration>(),
                   0.040816326530612244898 * SIUnit<Acceleration>(),
                   0.061224489795918367347 * SIUnit<Acceleration>()})));

  EXPECT_FALSE(fork->has_intrinsic_acceleration());
  fork->set_intrinsic_acceleration(
      [this](Instant const& t) {
        return Vector<Acceleration, World>(
            {4 * SIUnit<Length>() / ((t - t0_) * (t - t0_)),
             5 * SIUnit<Length>() / ((t - t0_) * (t - t0_)),
             6 * SIUnit<Length>() / ((t - t0_) * (t - t0_))});});
  EXPECT_TRUE(fork->has_intrinsic_acceleration());
  EXPECT_THAT(fork->evaluate_intrinsic_acceleration(t1_),
              Eq(Vector<Acceleration, World>({0 * SIUnit<Acceleration>(),
                                              0 * SIUnit<Acceleration>(),
                                              0 * SIUnit<Acceleration>()})));
  EXPECT_THAT(fork->evaluate_intrinsic_acceleration(t2_),
              Eq(Vector<Acceleration, World>({0 * SIUnit<Acceleration>(),
                                              0 * SIUnit<Acceleration>(),
                                              0 * SIUnit<Acceleration>()})));
  EXPECT_THAT(fork->evaluate_intrinsic_acceleration(t3_),
              Eq(Vector<Acceleration, World>(
                  {0.0054869684499314128944 * SIUnit<Acceleration>(),
                   0.0068587105624142661180 * SIUnit<Acceleration>(),
                   0.0082304526748971193416 * SIUnit<Acceleration>()})));
  fork->clear_intrinsic_acceleration();
  EXPECT_FALSE(fork->has_intrinsic_acceleration());
  EXPECT_TRUE(massless_trajectory_->has_intrinsic_acceleration());

  massless_trajectory_->clear_intrinsic_acceleration();
  EXPECT_FALSE(fork->has_intrinsic_acceleration());
  EXPECT_FALSE(massless_trajectory_->has_intrinsic_acceleration());
}

TEST_F(TrajectoryDeathTest, NativeIteratorError) {
  EXPECT_DEATH({
    Trajectory<World>::NativeIterator it = massive_trajectory_->last();
  }, "Empty trajectory");
  EXPECT_DEATH({
    Trajectory<World>::NativeIterator it = massive_trajectory_->first();
    ++it;
  }, "beyond end");
}

TEST_F(TrajectoryTest, NativeIteratorSuccess) {
  Trajectory<World>::NativeIterator it = massive_trajectory_->first();
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

  not_null<Trajectory<World>*> const fork = massless_trajectory_->Fork(t2_);
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

TEST_F(TrajectoryDeathTest, TransformingIteratorError) {
  EXPECT_DEATH({
    Trajectory<World>::TransformingIterator<World> it =
        massive_trajectory_->last_with_transform(massive_transform_);
  }, "Empty trajectory");
  EXPECT_DEATH({
    Trajectory<World>::TransformingIterator<World> it =
        massive_trajectory_->first_with_transform(massive_transform_);
    ++it;
  }, "beyond end");
}

TEST_F(TrajectoryTest, TransformingIteratorSuccess) {
  Trajectory<World>::TransformingIterator<World> it =
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

  not_null<Trajectory<World>*> const fork = massless_trajectory_->Fork(t2_);
  fork->Append(t4_, d4_);
  Trajectory<World>::Transform<World> const fork_transform =
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

TEST_F(TrajectoryTest, NativeIteratorOnOrAfterSuccess) {
  Trajectory<World>::NativeIterator it = massive_trajectory_->on_or_after(t0_);
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

  not_null<Trajectory<World>*> const fork = massless_trajectory_->Fork(t2_);
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

TEST_F(TrajectoryTest, TransformingIteratorOnOrAfterSuccess) {
  Trajectory<World>::TransformingIterator<World> it =
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

  not_null<Trajectory<World>*> const fork = massless_trajectory_->Fork(t2_);
  fork->Append(t4_, d4_);
  Trajectory<World>::Transform<World> const fork_transform =
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
