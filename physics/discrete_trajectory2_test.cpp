#include "physics/discrete_trajectory2.hpp"

#include <vector>

#include "base/serialization.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/serialization.hpp"
#include "testing_utilities/string_log_sink.hpp"

namespace principia {
namespace physics {

using base::ParseFromBytes;
using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Instant;
using geometry::Velocity;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::EqualsProto;
using testing_utilities::NewLinearTrajectoryTimeline;
using testing_utilities::ReadFromBinaryFile;
using testing_utilities::StringLogSink;
using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::HasSubstr;
using ::testing::Not;

class DiscreteTrajectory2Test : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;


  // Constructs a trajectory with three 5-second segments starting at |t0| and
  // the given |degrees_of_freedom|.
  DiscreteTrajectory2<World> MakeTrajectory(
      Instant const& t0,
      DegreesOfFreedom<World> const& degrees_of_freedom) {
    DiscreteTrajectory2<World> trajectory;
    std::optional<DegreesOfFreedom<World>> last_degrees_of_freedom;

    for (auto const& [t, degrees_of_freedom] :
         NewLinearTrajectoryTimeline(degrees_of_freedom,
                                     /*Δt=*/1 * Second,
                                     /*t1=*/t0,
                                     /*t2=*/t0 + 5 * Second)) {
      last_degrees_of_freedom = degrees_of_freedom;
      trajectory.Append(t, degrees_of_freedom);
    }

    trajectory.NewSegment();
    Velocity<World> const v2({0 * Metre / Second,
                              1 * Metre / Second,
                              0 * Metre / Second});
    for (auto const& [t, degrees_of_freedom] :
        NewLinearTrajectoryTimeline(DegreesOfFreedom<World>(
                                        last_degrees_of_freedom->position(),
                                        v2),
                                     /*Δt=*/1 * Second,
                                     /*t1=*/t0 + 5 * Second,
                                     /*t2=*/t0 + 10 * Second)) {
      last_degrees_of_freedom = degrees_of_freedom;
      trajectory.Append(t, degrees_of_freedom);
    }

    trajectory.NewSegment();
    Velocity<World> const v3({0 * Metre / Second,
                              0 * Metre / Second,
                              1 * Metre / Second});
    for (auto const& [t, degrees_of_freedom] :
        NewLinearTrajectoryTimeline(DegreesOfFreedom<World>(
                                        last_degrees_of_freedom->position(),
                                        v3),
                                     /*Δt=*/1 * Second,
                                     /*t1=*/t0 + 10 * Second,
                                     /*t2=*/t0 + 15 * Second)) {
      trajectory.Append(t, degrees_of_freedom);
    }

    return trajectory;
  }

  DiscreteTrajectory2<World> MakeTrajectory() {
    Velocity<World> const v1({1 * Metre / Second,
                              0 * Metre / Second,
                              0 * Metre / Second});
    return MakeTrajectory(t0_, DegreesOfFreedom<World>(World::origin, v1));
  }

  Instant const t0_;
};

TEST_F(DiscreteTrajectory2Test, Make) {
  auto const trajectory = MakeTrajectory();
}

TEST_F(DiscreteTrajectory2Test, IterateForward) {
  auto const trajectory = MakeTrajectory();
  std::vector<Instant> times;
  for (auto const& [t, _] : trajectory) {
    times.push_back(t);
  }
  EXPECT_THAT(times,
              ElementsAre(t0_,
                          t0_ + 1 * Second,
                          t0_ + 2 * Second,
                          t0_ + 3 * Second,
                          t0_ + 4 * Second,
                          t0_ + 5 * Second,
                          t0_ + 6 * Second,
                          t0_ + 7 * Second,
                          t0_ + 8 * Second,
                          t0_ + 9 * Second,
                          t0_ + 10 * Second,
                          t0_ + 11 * Second,
                          t0_ + 12 * Second,
                          t0_ + 13 * Second,
                          t0_ + 14 * Second));
}

TEST_F(DiscreteTrajectory2Test, IterateBackward) {
  auto const trajectory = MakeTrajectory();
  std::vector<Instant> times;
  for (auto it = trajectory.rbegin(); it != trajectory.rend(); ++it) {
    times.push_back(it->first);
  }
  EXPECT_THAT(times,
              ElementsAre(t0_ + 14 * Second,
                          t0_ + 13 * Second,
                          t0_ + 12 * Second,
                          t0_ + 11 * Second,
                          t0_ + 10 * Second,
                          t0_ + 9 * Second,
                          t0_ + 8 * Second,
                          t0_ + 7 * Second,
                          t0_ + 6 * Second,
                          t0_ + 5 * Second,
                          t0_ + 4 * Second,
                          t0_ + 3 * Second,
                          t0_ + 2 * Second,
                          t0_ + 1 * Second,
                          t0_));
}

TEST_F(DiscreteTrajectory2Test, Empty) {
  DiscreteTrajectory2<World> trajectory;
  EXPECT_TRUE(trajectory.empty());
  trajectory = MakeTrajectory();
  EXPECT_FALSE(trajectory.empty());
}

TEST_F(DiscreteTrajectory2Test, Size) {
  DiscreteTrajectory2<World> trajectory;
  EXPECT_EQ(0, trajectory.size());
  trajectory = MakeTrajectory();
  EXPECT_EQ(15, trajectory.size());
}

TEST_F(DiscreteTrajectory2Test, Find) {
  auto const trajectory = MakeTrajectory();
  {
    auto const it = trajectory.find(t0_ + 3 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 3 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({3 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory.find(t0_ + 13 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 13 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   4 * Metre,
                                                   3 * Metre}));
  }
  {
    auto const it = trajectory.find(t0_ + 3.14 * Second);
    EXPECT_TRUE(it == trajectory.end());
  }
}

TEST_F(DiscreteTrajectory2Test, LowerBound) {
  auto const trajectory = MakeTrajectory();
  {
    auto const it = trajectory.lower_bound(t0_ + 3.9 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 4 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory.lower_bound(t0_ + 4 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 4 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory.lower_bound(t0_ + 4.1 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 5 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory.lower_bound(t0_ + 13 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 13 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   4 * Metre,
                                                   3 * Metre}));
  }
  {
    auto const it = trajectory.lower_bound(t0_ + 14.2 * Second);
    EXPECT_TRUE(it == trajectory.end());
  }
}

TEST_F(DiscreteTrajectory2Test, UpperBound) {
  auto const trajectory = MakeTrajectory();
  {
    auto const it = trajectory.upper_bound(t0_ + 3.9 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 4 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory.upper_bound(t0_ + 4 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 5 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory.upper_bound(t0_ + 4.1 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 5 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory.upper_bound(t0_ + 13 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 14 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   4 * Metre,
                                                   4 * Metre}));
  }
  {
    auto const it = trajectory.upper_bound(t0_ + 14.2 * Second);
    EXPECT_TRUE(it == trajectory.end());
  }
}

TEST_F(DiscreteTrajectory2Test, Segments) {
  auto const trajectory = MakeTrajectory();
  std::vector<Instant> begin;
  std::vector<Instant> rbegin;
  for (auto const& sit : trajectory.segments()) {
    begin.push_back(sit.begin()->first);
    rbegin.push_back(sit.rbegin()->first);
  }
  EXPECT_THAT(
      begin,
      ElementsAre(t0_, t0_ + 4 * Second, t0_ + 9 * Second));
  EXPECT_THAT(
      rbegin,
      ElementsAre(t0_ + 4 * Second, t0_ + 9 * Second, t0_ + 14 * Second));
}

TEST_F(DiscreteTrajectory2Test, RSegments) {
  auto const trajectory = MakeTrajectory();
  std::vector<Instant> begin;
  std::vector<Instant> rbegin;
  for (auto const& sit : trajectory.rsegments()) {
    begin.push_back(sit.begin()->first);
    rbegin.push_back(sit.rbegin()->first);
  }
  EXPECT_THAT(
      begin,
      ElementsAre(t0_ + 9 * Second, t0_ + 4 * Second, t0_));
  EXPECT_THAT(
      rbegin,
      ElementsAre(t0_ + 14 * Second, t0_ + 9 * Second, t0_ + 4 * Second));
}

TEST_F(DiscreteTrajectory2Test, DetachSegments) {
  auto trajectory1 = MakeTrajectory();
  auto const first_segment = trajectory1.segments().begin();
  auto const second_segment = std::next(first_segment);
  auto trajectory2 = trajectory1.DetachSegments(second_segment);
  EXPECT_EQ(1, trajectory1.segments().size());
  EXPECT_EQ(2, trajectory2.segments().size());
  EXPECT_EQ(t0_, trajectory1.begin()->first);
  EXPECT_EQ(t0_ + 4 * Second, trajectory1.rbegin()->first);
  EXPECT_EQ(t0_ + 4 * Second, trajectory2.begin()->first);
  EXPECT_EQ(t0_ + 14 * Second, trajectory2.rbegin()->first);

  // Check that the trajectories are minimally usable (in particular, as far as
  // the time-to-segment mapping is concerned).
  {
    auto const it = trajectory1.lower_bound(t0_ + 3.9 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 4 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory1.lower_bound(t0_ + 4 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 4 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory2.lower_bound(t0_ + 4 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 4 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory2.lower_bound(t0_ + 4.1 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 5 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
}

TEST_F(DiscreteTrajectory2Test, AttachSegments) {
  auto trajectory1 = MakeTrajectory();
  auto trajectory2 = MakeTrajectory(
      t0_ + 14 * Second,
      DegreesOfFreedom<World>(
          World::origin + Displacement<World>({4 * Metre,
                                               4 * Metre,
                                               4 * Metre}),
          Velocity<World>({0 * Metre / Second,
                           0 * Metre / Second,
                           1 * Metre / Second})));
  trajectory1.AttachSegments(std::move(trajectory2));
  EXPECT_EQ(6, trajectory1.segments().size());
  EXPECT_EQ(t0_, trajectory1.begin()->first);
  EXPECT_EQ(t0_ + 28 * Second, trajectory1.rbegin()->first);

  // Check that the trajectories are minimally usable (in particular, as far as
  // the time-to-segment mapping is concerned).
  {
    auto const it = trajectory1.lower_bound(t0_ + 13.9 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 14 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   4 * Metre,
                                                   4 * Metre}));
  }
  {
    auto const it = trajectory1.lower_bound(t0_ + 14 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 14 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   4 * Metre,
                                                   4 * Metre}));
  }
  {
    auto const it = trajectory1.lower_bound(t0_ + 14.1 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 15 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   4 * Metre,
                                                   5 * Metre}));
  }
}

TEST_F(DiscreteTrajectory2Test, DeleteSegments) {
  auto trajectory = MakeTrajectory();
  auto const first_segment = trajectory.segments().begin();
  auto const second_segment = std::next(first_segment);
  trajectory.DeleteSegments(second_segment);
  EXPECT_EQ(1, trajectory.segments().size());
  EXPECT_EQ(t0_, trajectory.begin()->first);
  EXPECT_EQ(t0_ + 4 * Second, trajectory.rbegin()->first);
}

TEST_F(DiscreteTrajectory2Test, ForgetAfter) {
  auto trajectory = MakeTrajectory();

  trajectory.ForgetAfter(t0_ + 12 * Second);
  EXPECT_EQ(3, trajectory.segments().size());
  EXPECT_EQ(t0_, trajectory.begin()->first);
  EXPECT_EQ(t0_ + 11 * Second, trajectory.rbegin()->first);

  trajectory.ForgetAfter(t0_ + 6.1 * Second);
  EXPECT_EQ(2, trajectory.segments().size());
  EXPECT_EQ(t0_, trajectory.begin()->first);
  EXPECT_EQ(t0_ + 6 * Second, trajectory.rbegin()->first);

  trajectory.ForgetAfter(t0_ + 4 * Second);
  EXPECT_EQ(1, trajectory.segments().size());
  EXPECT_EQ(t0_, trajectory.begin()->first);
  EXPECT_EQ(t0_ + 4 * Second, trajectory.rbegin()->first);
}

TEST_F(DiscreteTrajectory2Test, ForgetBefore) {
  auto trajectory = MakeTrajectory();

  trajectory.ForgetBefore(t0_ + 3 * Second);
  EXPECT_EQ(3, trajectory.segments().size());
  EXPECT_EQ(t0_ + 3 * Second, trajectory.begin()->first);
  EXPECT_EQ(t0_ + 14 * Second, trajectory.rbegin()->first);

  trajectory.ForgetBefore(t0_ + 6.1 * Second);
  EXPECT_EQ(2, trajectory.segments().size());
  EXPECT_EQ(t0_ + 7 * Second, trajectory.begin()->first);
  EXPECT_EQ(t0_ + 14 * Second, trajectory.rbegin()->first);

  trajectory.ForgetBefore(t0_ + 9 * Second);
  EXPECT_EQ(1, trajectory.segments().size());
  EXPECT_EQ(t0_ + 9 * Second, trajectory.begin()->first);
  EXPECT_EQ(t0_ + 14 * Second, trajectory.rbegin()->first);
}

TEST_F(DiscreteTrajectory2Test, TMinTMaxEvaluate) {
  auto const trajectory = MakeTrajectory();
  EXPECT_EQ(t0_, trajectory.t_min());
  EXPECT_EQ(t0_ + 14 * Second, trajectory.t_max());
  EXPECT_THAT(trajectory.EvaluateDegreesOfFreedom(t0_ + 3.14 * Second),
      Componentwise(AlmostEquals(
                        World::origin + Displacement<World>({3.14 * Metre,
                                                             0 * Metre,
                                                             0 * Metre}), 0),
                    AlmostEquals(Velocity<World>({1 * Metre / Second,
                                                  0 * Metre / Second,
                                                  0 * Metre / Second}), 0)));
  EXPECT_THAT(trajectory.EvaluateDegreesOfFreedom(t0_ + 6.78 * Second),
      Componentwise(AlmostEquals(
                        World::origin + Displacement<World>({4 * Metre,
                                                             1.78 * Metre,
                                                             0 * Metre}), 1),
                    AlmostEquals(Velocity<World>({0 * Metre / Second,
                                        1 * Metre / Second,
                                        0 * Metre / Second}), 0)));
}

TEST_F(DiscreteTrajectory2Test, SerializationRoundTrip) {
  auto const trajectory = MakeTrajectory();
  auto const trajectory_first_segment = trajectory.segments().begin();
  auto const trajectory_second_segment = std::next(trajectory_first_segment);

  serialization::DiscreteTrajectory message1;
  trajectory.WriteToMessage(&message1,
                            /*tracked=*/{trajectory_second_segment},
                            /*exact=*/
                            {trajectory.lower_bound(t0_ + 2 * Second),
                             trajectory.lower_bound(t0_ + 3 * Second)});

  DiscreteTrajectorySegmentIterator<World> deserialized_second_segment;
  auto const deserialized_trajectory =
      DiscreteTrajectory2<World>::ReadFromMessage(
          message1, /*tracked=*/{&deserialized_second_segment});

  // Check that the tracked segment was properly retrieved.
  EXPECT_EQ(t0_ + 4 * Second, deserialized_second_segment->begin()->first);
  EXPECT_EQ(t0_ + 9 * Second, deserialized_second_segment->rbegin()->first);

  // Check that the exact points are exact.
  EXPECT_EQ(deserialized_trajectory.lower_bound(t0_ + 2 * Second)->second,
            trajectory.lower_bound(t0_ + 2 * Second)->second);
  EXPECT_EQ(deserialized_trajectory.lower_bound(t0_ + 3 * Second)->second,
            trajectory.lower_bound(t0_ + 3 * Second)->second);

  serialization::DiscreteTrajectory message2;
  deserialized_trajectory.WriteToMessage(
      &message2,
      /*tracked=*/{deserialized_second_segment},
      /*exact=*/
      {deserialized_trajectory.lower_bound(t0_ + 2 * Second),
       deserialized_trajectory.lower_bound(t0_ + 3 * Second)});

  EXPECT_THAT(message2, EqualsProto(message1));
}

TEST_F(DiscreteTrajectory2Test, SerializationPreHaarCompatibility) {
  StringLogSink log_warning(google::WARNING);
  auto const serialized_message = ReadFromBinaryFile(
      R"(P:\Public Mockingbird\Principia\Saves\3136\trajectory_3136.proto.bin)");  // NOLINT
  auto const message =
      ParseFromBytes<serialization::DiscreteTrajectory>(serialized_message);
  DiscreteTrajectory2<World>::SegmentIterator psychohistory;
  auto const history = DiscreteTrajectory2<World>::ReadFromMessage(
      message, /*tracked=*/{&psychohistory});
  EXPECT_THAT(log_warning.string(),
              AllOf(HasSubstr("pre-Ζήνων"), Not(HasSubstr("pre-Haar"))));

  // Note that the sizes don't have the same semantics as pre-Haar.
  EXPECT_EQ(435'929, history.size());
  EXPECT_EQ(2, psychohistory->size());
}

}  // namespace physics
}  // namespace principia
