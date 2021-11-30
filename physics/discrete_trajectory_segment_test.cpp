#include "physics/discrete_trajectory_segment.hpp"

#include <algorithm>
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
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace physics {

using base::make_not_null_unique;
using base::not_null;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Instant;
using geometry::Velocity;
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
using testing_utilities::AbsoluteErrorFrom;
using testing_utilities::AlmostEquals;
using testing_utilities::AppendTrajectoryTimeline;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::NewCircularTrajectoryTimeline;
using testing_utilities::NewLinearTrajectoryTimeline;
using testing_utilities::operator""_⑴;
using ::testing::Eq;
using ::testing::Le;
using ::testing::Lt;

class DiscreteTrajectorySegmentTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using Segments = internal_discrete_trajectory_types::Segments<World>;

  DiscreteTrajectorySegmentTest()
      : segments_(MakeSegments(1)) {
    segment_ = &*segments_->begin();

    EXPECT_OK(segment_->Append(t0_ + 2 * Second, unmoving_origin_));
    EXPECT_OK(segment_->Append(t0_ + 3 * Second, unmoving_origin_));
    EXPECT_OK(segment_->Append(t0_ + 5 * Second, unmoving_origin_));
    EXPECT_OK(segment_->Append(t0_ + 7 * Second, unmoving_origin_));
    EXPECT_OK(segment_->Append(t0_ + 11 * Second, unmoving_origin_));
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

  void ForgetBefore(Instant const& t,
                    DiscreteTrajectorySegment<World>& segment) {
    segment.ForgetBefore(t);
  }

  static DiscreteTrajectorySegmentIterator<World> MakeIterator(
      not_null<Segments*> const segments,
      typename Segments::iterator iterator) {
    return DiscreteTrajectorySegmentIterator<World>(segments, iterator);
  }

  // Constructs a list of |n| segments which are properly initialized.
  // TODO(phl): Move to a central place.
  static not_null<std::unique_ptr<Segments>> MakeSegments(const int n) {
    auto segments = make_not_null_unique<Segments>(n);
    for (auto it = segments->begin(); it != segments->end(); ++it) {
      *it = DiscreteTrajectorySegment<World>(MakeIterator(segments.get(), it));
    }
    return segments;
  }

  DiscreteTrajectorySegment<World>* segment_;
  not_null<std::unique_ptr<Segments>> segments_;
  Instant const t0_;
  DegreesOfFreedom<World> unmoving_origin_{World::origin, World::unmoving};
};

TEST_F(DiscreteTrajectorySegmentTest, BackFront) {
  EXPECT_EQ(t0_ + 2 * Second, segment_->front().time);
  EXPECT_EQ(t0_ + 11 * Second, segment_->back().time);
}

TEST_F(DiscreteTrajectorySegmentTest, Extremities) {
  {
    auto const it = segment_->begin();
    EXPECT_EQ(t0_ + 2 * Second, it->time);
  }
  {
    auto it = segment_->end();
    --it;
    EXPECT_EQ(t0_ + 11 * Second, it->time);
  }
  {
    auto const it = segment_->rbegin();
    EXPECT_EQ(t0_ + 11 * Second, it->time);
  }
  {
    auto it = segment_->rend();
    --it;
    EXPECT_EQ(t0_ + 2 * Second, it->time);
  }
}

TEST_F(DiscreteTrajectorySegmentTest, Find) {
  {
    auto const it = segment_->find(t0_ + 5 * Second);
    EXPECT_EQ(t0_ + 5 * Second, it->time);
  }
  {
    auto const it = segment_->find(t0_ + 6 * Second);
    EXPECT_TRUE(it == segment_->end());
  }
}

TEST_F(DiscreteTrajectorySegmentTest, LowerBoundUpperBound) {
  {
    auto const it = segment_->lower_bound(t0_ + 5 * Second);
    EXPECT_EQ(t0_ + 5 * Second, it->time);
  }
  {
    auto const it = segment_->lower_bound(t0_ + 6 * Second);
    EXPECT_EQ(t0_ + 7 * Second, it->time);
  }
  {
    auto const it = segment_->lower_bound(t0_ + 12 * Second);
    EXPECT_TRUE(it == segment_->end());
  }
  {
    auto const it = segment_->upper_bound(t0_ + 5 * Second);
    EXPECT_EQ(t0_ + 7 * Second, it->time);
  }
  {
    auto const it = segment_->upper_bound(t0_ + 6 * Second);
    EXPECT_EQ(t0_ + 7 * Second, it->time);
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
  EXPECT_EQ(t0_ + 3 * Second, segment_->rbegin()->time);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetAfterNonexisting) {
  ForgetAfter(t0_ + 6 * Second);
  EXPECT_EQ(t0_ + 5 * Second, segment_->rbegin()->time);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetAfterPastTheEnd) {
  ForgetAfter(t0_ + 29 * Second);
  EXPECT_EQ(t0_ + 11 * Second, segment_->rbegin()->time);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetBeforeExisting) {
  ForgetBefore(t0_ + 7 * Second);
  EXPECT_EQ(t0_ + 7 * Second, segment_->begin()->time);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetBeforeNonexisting) {
  ForgetBefore(t0_ + 6 * Second);
  EXPECT_EQ(t0_ + 7 * Second, segment_->begin()->time);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetBeforeTheBeginning) {
  ForgetBefore(t0_ + 1 * Second);
  EXPECT_EQ(t0_ + 2 * Second, segment_->begin()->time);
}

TEST_F(DiscreteTrajectorySegmentTest, Evaluate) {
  auto const segments = MakeSegments(1);
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

TEST_F(DiscreteTrajectorySegmentTest, DownsamplingCircle) {
  auto const circle_segments = MakeSegments(1);
  auto const downsampled_circle_segments = MakeSegments(1);
  auto& circle = *circle_segments->begin();
  auto& downsampled_circle = *downsampled_circle_segments->begin();
  downsampled_circle.SetDownsampling(
      {.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)});
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

TEST_F(DiscreteTrajectorySegmentTest, DownsamplingStraightLine) {
  auto const line_segments = MakeSegments(1);
  auto const downsampled_line_segments = MakeSegments(1);
  auto& line = *line_segments->begin();
  auto& downsampled_line = *downsampled_line_segments->begin();
  downsampled_line.SetDownsampling(
      {.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)});
  auto const v = Velocity<World>({1 * Metre / Second,
                                  2 * Metre / Second,
                                  3 * Metre / Second});
  Time const Δt = 10 * Milli(Second);
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 10 * Second;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline<World>(v, Δt, t1, t2),
      /*to=*/line);
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline<World>(v, Δt, t1, t2),
      /*to=*/downsampled_line);

  EXPECT_THAT(line.size(), Eq(1001));
  // In the test3200 release this used to have 1001 points, see #3203.
  EXPECT_THAT(downsampled_line.size(), Eq(21));
  std::vector<Length> position_errors;
  std::vector<Speed> velocity_errors;
  for (auto const& [time, degrees_of_freedom] : line) {
    position_errors.push_back(
        (downsampled_line.EvaluatePosition(time) -
         degrees_of_freedom.position()).Norm());
    velocity_errors.push_back(
        (downsampled_line.EvaluateVelocity(time) -
         degrees_of_freedom.velocity()).Norm());
  }
  EXPECT_THAT(*std::max_element(position_errors.begin(), position_errors.end()),
              IsNear(3.6e-15_⑴ * Metre));
  EXPECT_THAT(*std::max_element(velocity_errors.begin(), velocity_errors.end()),
              IsNear(1.1e-14_⑴ * Metre / Second));
}

TEST_F(DiscreteTrajectorySegmentTest, DownsamplingForgetAfter) {
  auto const circle_segments = MakeSegments(1);
  auto const forgotten_circle_segments = MakeSegments(1);
  auto& circle = *circle_segments->begin();
  auto& forgotten_circle = *forgotten_circle_segments->begin();
  circle.SetDownsampling(
      {.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)});
  forgotten_circle.SetDownsampling(
      {.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)});
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
  Instant const restart_time = forgotten_circle.lower_bound(t2)->time;
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
        (circle.find(t)->degrees_of_freedom.position() -
         degrees_of_freedom.position()).Norm());
    velocity_errors.push_back(
        (circle.find(t)->degrees_of_freedom.velocity() -
         degrees_of_freedom.velocity()).Norm());
  }
  EXPECT_THAT(*std::max_element(position_errors.begin(), position_errors.end()),
              AlmostEquals(0 * Metre, 0));
  EXPECT_THAT(*std::max_element(velocity_errors.begin(), velocity_errors.end()),
              AlmostEquals(0 * Metre / Second, 0));
}

TEST_F(DiscreteTrajectorySegmentTest, SerializationWithDownsampling) {
  auto const circle_segments = MakeSegments(1);
  auto& circle = *circle_segments->begin();
  circle.SetDownsampling(
      {.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)});
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 10 * Milli(Second);
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 5 * Second;
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, t1, t2),
      /*to=*/circle);
  auto const circle_t_max = circle.t_max();

  serialization::DiscreteTrajectorySegment message;
  circle.WriteToMessage(&message,
                        /*exact=*/{circle.lower_bound(t0_ + 2 * Second),
                                   circle.lower_bound(t0_ + 3 * Second)});

  auto const deserialized_circle_segments = MakeSegments(1);
  auto& deserialized_circle = *deserialized_circle_segments->begin();
  deserialized_circle =
      DiscreteTrajectorySegment<World>::ReadFromMessage(
          message,
          /*self=*/MakeIterator(deserialized_circle_segments.get(),
                                deserialized_circle_segments->begin()));

  // Serialization/deserialization preserves the size, the times, and nudges the
  // positions by less than the tolerance.  It also preserve the degrees of
  // freedom at the "exact" iterators.
  EXPECT_THAT(circle.size(), Eq(39));
  EXPECT_THAT(deserialized_circle.size(), circle.size());
  for (auto it1 = circle.begin(), it2 = deserialized_circle.begin();
       it1 != circle.end();
       ++it1, ++it2) {
    auto const& [t1, degrees_of_freedom1] = *it1;
    auto const& [t2, degrees_of_freedom2] = *it2;
    EXPECT_EQ(t1, t2);
    EXPECT_THAT(degrees_of_freedom2.position(),
                AbsoluteErrorFrom(degrees_of_freedom1.position(),
                                  Lt(0.22 * Milli(Metre))));
    EXPECT_THAT(degrees_of_freedom2.velocity(),
                AbsoluteErrorFrom(degrees_of_freedom1.velocity(),
                                  Lt(1.1 * Milli(Metre) / Second)));
  }
  EXPECT_NE(
      deserialized_circle.lower_bound(t0_ + 1 * Second)->degrees_of_freedom,
      circle.lower_bound(t0_ + 1 * Second)->degrees_of_freedom);
  EXPECT_EQ(
      deserialized_circle.lower_bound(t0_ + 2 * Second)->degrees_of_freedom,
      circle.lower_bound(t0_ + 2 * Second)->degrees_of_freedom);
  EXPECT_EQ(
      deserialized_circle.lower_bound(t0_ + 3 * Second)->degrees_of_freedom,
      circle.lower_bound(t0_ + 3 * Second)->degrees_of_freedom);
  EXPECT_NE(
      deserialized_circle.lower_bound(t0_ + 4 * Second)->degrees_of_freedom,
      circle.lower_bound(t0_ + 4 * Second)->degrees_of_freedom);

  // Appending may result in different downsampling because the positions differ
  // a bit.
  Instant const t3 = t0_ + 10 * Second;
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, circle_t_max + Δt, t3),
      /*to=*/circle);
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, circle_t_max + Δt, t3),
      /*to=*/deserialized_circle);

  // Despite the difference in downsampling (and therefore in size) the two
  // trajectories are still within the tolerance.
  EXPECT_THAT(circle.size(), Eq(77));
  EXPECT_THAT(deserialized_circle.size(), Eq(78));
  for (Instant t = t0_;
       t <= std::min(circle.rbegin()->time, deserialized_circle.rbegin()->time);
       t += Δt) {
    EXPECT_THAT(
        deserialized_circle.EvaluatePosition(t),
        AbsoluteErrorFrom(circle.EvaluatePosition(t), Le(0.23 * Milli(Metre))));
    EXPECT_THAT(deserialized_circle.EvaluateVelocity(t),
                AbsoluteErrorFrom(circle.EvaluateVelocity(t),
                                  Le(5.7 * Milli(Metre) / Second)));
  }
}

TEST_F(DiscreteTrajectorySegmentTest, SerializationRoundTrip) {
  auto const circle_segments = MakeSegments(1);
  auto& circle = *circle_segments->begin();
  circle.SetDownsampling(
      {.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)});
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 10 * Milli(Second);
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 5 * Second;
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, t1, t2),
      /*to=*/circle);

  serialization::DiscreteTrajectorySegment message1;
  circle.WriteToMessage(
      &message1,
      /*exact=*/{circle.lower_bound(t0_ + 2 * Second),
                 circle.lower_bound(t0_ + 3 * Second)});

  auto const deserialized_circle_segments = MakeSegments(1);
  auto& deserialized_circle = *deserialized_circle_segments->begin();
  deserialized_circle =
      DiscreteTrajectorySegment<World>::ReadFromMessage(
          message1,
          /*self=*/MakeIterator(deserialized_circle_segments.get(),
                                deserialized_circle_segments->begin()));

  serialization::DiscreteTrajectorySegment message2;
  deserialized_circle.WriteToMessage(
      &message2,
      /*exact=*/{circle.lower_bound(t0_ + 2 * Second),
                 circle.lower_bound(t0_ + 3 * Second)});

  EXPECT_THAT(message2, EqualsProto(message1));
}

TEST_F(DiscreteTrajectorySegmentTest, SerializationEmpty) {
  DiscreteTrajectorySegment<World> segment;
  serialization::DiscreteTrajectorySegment message;
  segment.WriteToMessage(&message, /*exact=*/{});
  auto const deserialized_segments = MakeSegments(1);
  auto& deserialized_segment = *deserialized_segments->begin();
  deserialized_segment = DiscreteTrajectorySegment<World>::ReadFromMessage(
      message,
      /*self=*/MakeIterator(deserialized_segments.get(),
                            deserialized_segments->begin()));
}

TEST_F(DiscreteTrajectorySegmentTest, SerializationRange) {
  auto const circle1_segments = MakeSegments(1);
  auto const circle2_segments = MakeSegments(1);
  auto& circle1 = *circle1_segments->begin();
  auto& circle2 = *circle2_segments->begin();
  circle1.SetDownsampling(
      {.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)});
  circle2.SetDownsampling(
      {.max_dense_intervals = 50, .tolerance = 1 * Milli(Metre)});
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 10 * Milli(Second);
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 5 * Second;
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, t1, t2),
      /*to=*/circle1);
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(ω, r, Δt, t1, t2),
      /*to=*/circle2);

  serialization::DiscreteTrajectorySegment message1;
  circle1.WriteToMessage(&message1,
                         /*begin=*/circle1.upper_bound(t0_ + 4.9 * Second),
                         /*end=*/circle1.upper_bound(t0_ + 4.98 * Second),
                         /*exact=*/{});
  EXPECT_LE(message1.number_of_dense_points(), message1.zfp().timeline_size());

  serialization::DiscreteTrajectorySegment message2;
  ForgetBefore(circle2.upper_bound(t0_ + 4.9 * Second)->time, circle2);
  ForgetAfter(circle2.upper_bound(t0_ + 4.98 * Second)->time, circle2);
  circle2.WriteToMessage(&message2,
                         /*exact=*/{});
  EXPECT_LE(message2.number_of_dense_points(), message2.zfp().timeline_size());

  // Writing a range of the segment is equivalent to forgetting and writing the
  // result.
  EXPECT_THAT(message1, EqualsProto(message2));
}

}  // namespace physics
}  // namespace principia
