#pragma once

#include <cstdint>
#include <iterator>
#include <optional>

#include "absl/container/btree_map.h"
#include "absl/status/status.h"
#include "geometry/named_quantities.hpp"
#include "numerics/hermite3.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/trajectory.hpp"

namespace principia {

namespace testing_utilities {
FORWARD_DECLARE_FROM(discrete_trajectory_factories,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectoryFactoriesFriend);
}  // namespace testing_utilities

namespace physics {

FORWARD_DECLARE_FROM(discrete_trajectory2,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectory2);

class DiscreteTrajectoryIteratorTest;
class DiscreteTrajectorySegmentIteratorTest;
class DiscreteTrajectorySegmentTest;

namespace internal_discrete_trajectory_segment {

using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using numerics::Hermite3;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTrajectorySegment : public Trajectory<Frame> {
  using Timeline = internal_discrete_trajectory_types::Timeline<Frame>;

 public:
  using key_type = typename Timeline::key_type;
  using value_type = typename Timeline::value_type;

  using iterator = DiscreteTrajectoryIterator<Frame>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  // TODO(phl): Decide which constructors should be public.
  DiscreteTrajectorySegment() = default;
  explicit DiscreteTrajectorySegment(
      DiscreteTrajectorySegmentIterator<Frame> self);

  ~DiscreteTrajectorySegment() = default;

  // Moveable.
  DiscreteTrajectorySegment(DiscreteTrajectorySegment&&) = default;
  DiscreteTrajectorySegment& operator=(DiscreteTrajectorySegment&&) = default;
  DiscreteTrajectorySegment(const DiscreteTrajectorySegment&) = delete;
  DiscreteTrajectorySegment& operator=(const DiscreteTrajectorySegment&) =
      delete;

  iterator begin() const;
  iterator end() const;

  reverse_iterator rbegin() const;
  reverse_iterator rend() const;

  // TODO(phl): We probably don't want empty segments.
  bool empty() const;
  std::int64_t size() const;

  iterator find(Instant const& t) const;

  iterator lower_bound(Instant const& t) const;
  iterator upper_bound(Instant const& t) const;

  Instant t_min() const override;
  Instant t_max() const override;

  Position<Frame> EvaluatePosition(Instant const& t) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& t) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& t) const override;

 private:
  using DownsamplingParameters =
      internal_discrete_trajectory_types::DownsamplingParameters;

  // Changes the |self_| iterator.  Only for use when attaching/detaching
  // segments.
  void SetSelf(DiscreteTrajectorySegmentIterator<Frame> self);

  absl::Status Append(Instant const& t,
                      DegreesOfFreedom<Frame> const& degrees_of_freedom);

  // Removes all points with a time greater than or equal to |t| (1st overload)
  // or starting at |begin| (2nd overload).
  void ForgetAfter(Instant const& t);
  void ForgetAfter(typename Timeline::const_iterator begin);

  // Removes all points with a time strictly less than |t| (1st overload) or
  // ending at |end| (2nd overload).
  void ForgetBefore(Instant const& t);
  void ForgetBefore(typename Timeline::const_iterator end);

  // This segment must have 0 or 1 points.  Occasionally removes intermediate
  // points from the segment when |Append|ing, ensuring that positions remain
  // within the desired tolerance.
  void SetDownsampling(DownsamplingParameters const& downsampling_parameters);

  // Clear the downsampling parameters.  From now on, all points appended to the
  // segment are going to be retained.
  void ClearDownsampling();

  // Called by |Append| after appending a point to this segment.  If
  // appropriate, performs downsampling and deletes some of the points of the
  // segment.
  absl::Status DownsampleIfNeeded();

  // Returns the Hermite interpolation for the left-open, right-closed
  // trajectory segment bounded above by |upper|.
  Hermite3<Instant, Position<Frame>> GetInterpolation(
      typename Timeline::const_iterator upper) const;

  typename Timeline::const_iterator timeline_begin() const;
  typename Timeline::const_iterator timeline_end() const;
  bool timeline_empty() const;

  std::optional<DownsamplingParameters> downsampling_parameters_;

  DiscreteTrajectorySegmentIterator<Frame> self_;
  Timeline timeline_;

  // The number of points at the end of the segment that are part of a "dense"
  // span, i.e., have not been subjected to downsampling yet.
  std::int64_t number_of_dense_points_ = 0;

  template<typename F>
  friend class physics::DiscreteTrajectory2;
  template<typename F>
  friend class physics::DiscreteTrajectoryIterator;

  // For testing.
  friend class physics::DiscreteTrajectoryIteratorTest;
  friend class physics::DiscreteTrajectorySegmentIteratorTest;
  friend class physics::DiscreteTrajectorySegmentTest;
  template<typename F>
  friend class testing_utilities::DiscreteTrajectoryFactoriesFriend;
};

}  // namespace internal_discrete_trajectory_segment

using internal_discrete_trajectory_segment::DiscreteTrajectorySegment;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_segment_body.hpp"
