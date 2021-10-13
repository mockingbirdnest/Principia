#include "physics/discrete_trajectory_segment.hpp"

#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/discrete_trajectory_types.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace physics {

using base::make_not_null_unique;
using base::not_null;
using geometry::Frame;
using geometry::Instant;
using quantities::Abs;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Speed;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Nano;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::AppendTrajectoryTimeline;
using testing_utilities::IsNear;
using testing_utilities::NewCircularTrajectoryTimeline;
using testing_utilities::operator""_⑴;
using ::testing::Eq;

class DiscreteTrajectorySegmentTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  using Segments = internal_discrete_trajectory_types::Segments<World>;

  DiscreteTrajectorySegmentTest()
      : segments_(MakeSegments()) {
    segment_ = &*segments_->begin();

    segment_->Append(t0_ + 2 * Second, unmoving_origin_);
    segment_->Append(t0_ + 3 * Second, unmoving_origin_);
    segment_->Append(t0_ + 5 * Second, unmoving_origin_);
    segment_->Append(t0_ + 7 * Second, unmoving_origin_);
    segment_->Append(t0_ + 11 * Second, unmoving_origin_);
  }

  // Constructs a list of one segment which is properly initialized.
  static not_null<std::unique_ptr<Segments>> MakeSegments() {
    auto segments = make_not_null_unique<Segments>(1);
    auto it = segments->begin();
    *it = DiscreteTrajectorySegment<World>(
        DiscreteTrajectorySegmentIterator<World>(segments.get(), it));
    return segments;
  }

  void ForgetAfter(Instant const& t) {
    segment_->ForgetAfter(t);
  }

  void ForgetAfter(Instant const& t,
                   DiscreteTrajectorySegment<World>& segment) {
    segment.ForgetAfter(t);
  }

  void ForgetBefore(Instant const& t) {
    segment_->ForgetBefore(t);
  }

  void SetDownsampling(
      internal_discrete_trajectory_types::DownsamplingParameters const&
          downsampling_parameters,
      DiscreteTrajectorySegment<World>& segment) {
    segment.SetDownsampling(downsampling_parameters);
  }

  DiscreteTrajectorySegment<World>* segment_;
  not_null<std::unique_ptr<Segments>> segments_;
  Instant const t0_;
  DegreesOfFreedom<World> unmoving_origin_{World::origin, World::unmoving};
};

TEST_F(DiscreteTrajectorySegmentTest, Extremities) {
  {
    auto const it = segment_->begin();
    EXPECT_EQ(t0_ + 2 * Second, it->first);
  }
  {
    auto it = segment_->end();
    --it;
    EXPECT_EQ(t0_ + 11 * Second, it->first);
  }
  {
    auto const it = segment_->rbegin();
    EXPECT_EQ(t0_ + 11 * Second, it->first);
  }
  {
    auto it = segment_->rend();
    --it;
    EXPECT_EQ(t0_ + 2 * Second, it->first);
  }
}

TEST_F(DiscreteTrajectorySegmentTest, Find) {
  {
    auto const it = segment_->find(t0_ + 5 * Second);
    EXPECT_EQ(t0_ + 5 * Second, it->first);
  }
  {
    auto const it = segment_->find(t0_ + 6 * Second);
    EXPECT_TRUE(it == segment_->end());
  }
}

TEST_F(DiscreteTrajectorySegmentTest, LowerBoundUpperBound) {
  {
    auto const it = segment_->lower_bound(t0_ + 5 * Second);
    EXPECT_EQ(t0_ + 5 * Second, it->first);
  }
  {
    auto const it = segment_->lower_bound(t0_ + 6 * Second);
    EXPECT_EQ(t0_ + 7 * Second, it->first);
  }
  {
    auto const it = segment_->lower_bound(t0_ + 12 * Second);
    EXPECT_TRUE(it == segment_->end());
  }
  {
    auto const it = segment_->upper_bound(t0_ + 5 * Second);
    EXPECT_EQ(t0_ + 7 * Second, it->first);
  }
  {
    auto const it = segment_->upper_bound(t0_ + 6 * Second);
    EXPECT_EQ(t0_ + 7 * Second, it->first);
  }
  {
    auto const it = segment_->upper_bound(t0_ + 11 * Second);
    EXPECT_TRUE(it == segment_->end());
  }
}

TEST_F(DiscreteTrajectorySegmentTest, EmptySize) {
  EXPECT_FALSE(segment_->empty());
  EXPECT_EQ(5, segment_->size());
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetAfterExisting) {
  ForgetAfter(t0_ + 5 * Second);
  EXPECT_EQ(t0_ + 3 * Second, segment_->rbegin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetAfterNonexisting) {
  ForgetAfter(t0_ + 6 * Second);
  EXPECT_EQ(t0_ + 5 * Second, segment_->rbegin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetAfterPastTheEnd) {
  ForgetAfter(t0_ + 29 * Second);
  EXPECT_EQ(t0_ + 11 * Second, segment_->rbegin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetBeforeExisting) {
  ForgetBefore(t0_ + 7 * Second);
  EXPECT_EQ(t0_ + 7 * Second, segment_->begin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetBeforeNonexisting) {
  ForgetBefore(t0_ + 6 * Second);
  EXPECT_EQ(t0_ + 7 * Second, segment_->begin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetBeforeTheBeginning) {
  ForgetBefore(t0_ + 1 * Second);
  EXPECT_EQ(t0_ + 2 * Second, segment_->begin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, Evaluate) {
  auto const segments = MakeSegments();
  auto& circle = *segments->begin();
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 10 * Milli(Second);
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 10 * Second;
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, t1, t2), /*to=*/circle);

  EXPECT_THAT(circle.size(), Eq(1001));
  std::vector<Length> position_errors;
  std::vector<Speed> velocity_errors;
  for (Instant t = circle.t_min();
       t <= circle.t_max();
       t += 1 * Milli(Second)) {
    position_errors.push_back(
        Abs((circle.EvaluatePosition(t) - World::origin).Norm() - r));
    velocity_errors.push_back(
        Abs(circle.EvaluateVelocity(t).Norm() - r * ω / Radian));
  }
  EXPECT_THAT(*std::max_element(position_errors.begin(), position_errors.end()),
              IsNear(4.2_⑴ * Nano(Metre)));
  EXPECT_THAT(*std::max_element(velocity_errors.begin(), velocity_errors.end()),
              IsNear(10.4_⑴ * Nano(Metre / Second)));
}

TEST_F(DiscreteTrajectorySegmentTest, Downsampling) {
  auto const circle_segments = MakeSegments();
  auto const downsampled_circle_segments = MakeSegments();
  auto& circle = *circle_segments->begin();
  auto& downsampled_circle = *downsampled_circle_segments->begin();
  SetDownsampling({.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)},
                  downsampled_circle);
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 10 * Milli(Second);
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 10 * Second;
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, t1, t2),
      /*to=*/circle);
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, t1, t2),
      /*to=*/downsampled_circle);

  EXPECT_THAT(circle.size(), Eq(1001));
  EXPECT_THAT(downsampled_circle.size(), Eq(77));
  std::vector<Length> position_errors;
  std::vector<Speed> velocity_errors;
  for (auto const& [time, degrees_of_freedom] : circle) {
    position_errors.push_back(
        (downsampled_circle.EvaluatePosition(time) -
         degrees_of_freedom.position()).Norm());
    velocity_errors.push_back(
        (downsampled_circle.EvaluateVelocity(time) -
         degrees_of_freedom.velocity()).Norm());
  }
  EXPECT_THAT(*std::max_element(position_errors.begin(), position_errors.end()),
              IsNear(0.98_⑴ * Milli(Metre)));
  EXPECT_THAT(*std::max_element(velocity_errors.begin(), velocity_errors.end()),
              IsNear(14_⑴ * Milli(Metre / Second)));
}

TEST_F(DiscreteTrajectorySegmentTest, DownsamplingForgetAfter) {
  auto const circle_segments = MakeSegments();
  auto const forgotten_circle_segments = MakeSegments();
  auto& circle = *circle_segments->begin();
  auto& forgotten_circle = *forgotten_circle_segments->begin();
  SetDownsampling({.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)},
                  circle);
  SetDownsampling({.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)},
                  forgotten_circle);
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 1.0 / 128.0 * Second;  // Yields exact times.
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 5 * Second;
  Instant const t3 = t0_ + 10 * Second;

  // Construct two identical trajectories with downsampling.
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, t1, t3),
      /*to=*/circle);
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, t1, t3),
      /*to=*/forgotten_circle);

  // Forget one of the trajectories in the middle, and append new points.
  Instant const restart_time = forgotten_circle.lower_bound(t2)->first;
  ForgetAfter(t2, forgotten_circle);
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, restart_time, t3),
      /*to=*/forgotten_circle);

  EXPECT_THAT(circle.size(), Eq(92));
  EXPECT_THAT(forgotten_circle.size(), Eq(circle.size()));
  std::vector<Length> position_errors;
  std::vector<Speed> velocity_errors;

  // Check that the two trajectories are identical.
  for (auto const [t, degrees_of_freedom] : forgotten_circle) {
    position_errors.push_back(
        (circle.find(t)->second.position() -
         degrees_of_freedom.position()).Norm());
    velocity_errors.push_back(
        (circle.find(t)->second.velocity() -
         degrees_of_freedom.velocity()).Norm());
  }
  EXPECT_THAT(*std::max_element(position_errors.begin(), position_errors.end()),
              AlmostEquals(0 * Metre, 0));
  EXPECT_THAT(*std::max_element(velocity_errors.begin(), velocity_errors.end()),
              AlmostEquals(0 * Metre / Second, 0));
}

}  // namespace physics
}  // namespace principia
