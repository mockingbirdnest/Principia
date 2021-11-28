#pragma once

#include <iterator>
#include <list>
#include <memory>
#include <vector>

#include "absl/container/btree_map.h"
#include "absl/status/status.h"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "base/tags.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_segment_range.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/trajectory.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {

FORWARD_DECLARE_FROM(discrete_trajectory_segment,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectorySegment);

namespace internal_discrete_trajectory {

using base::not_null;
using base::uninitialized_t;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTrajectory : public Trajectory<Frame> {
 public:
  using key_type =
      typename internal_discrete_trajectory_types::Timeline<Frame>::key_type;
  using value_type =
      typename internal_discrete_trajectory_types::Timeline<Frame>::value_type;

  using iterator = DiscreteTrajectoryIterator<Frame>;
  using reference = value_type const&;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using SegmentIterator = DiscreteTrajectorySegmentIterator<Frame>;
  using ReverseSegmentIterator = std::reverse_iterator<SegmentIterator>;
  using SegmentRange = DiscreteTrajectorySegmentRange<SegmentIterator>;
  using ReverseSegmentRange =
      DiscreteTrajectorySegmentRange<ReverseSegmentIterator>;

  DiscreteTrajectory();

  // Moveable.
  DiscreteTrajectory(DiscreteTrajectory&&) = default;
  DiscreteTrajectory& operator=(DiscreteTrajectory&&) = default;
  DiscreteTrajectory(const DiscreteTrajectory&) = delete;
  DiscreteTrajectory& operator=(const DiscreteTrajectory&) = delete;

  reference front() const;
  reference back() const;

  iterator begin() const;
  iterator end() const;

  reverse_iterator rbegin() const;
  reverse_iterator rend() const;

  bool empty() const;
  std::int64_t size() const;

  // Doesn't invalidate iterators to the first segment.
  void clear();

  iterator find(Instant const& t) const;

  iterator lower_bound(Instant const& t) const;
  iterator upper_bound(Instant const& t) const;

  SegmentRange segments() const;
  // TODO(phl): In C++20 this should be a reverse_view on segments.
  ReverseSegmentRange rsegments() const;

  SegmentIterator NewSegment();

  DiscreteTrajectory DetachSegments(SegmentIterator begin);
  SegmentIterator AttachSegments(DiscreteTrajectory trajectory);
  void DeleteSegments(SegmentIterator& begin);

  // Deletes the trajectory points with a time in [t, end[.  Drops the segments
  // that are empty as a result.
  void ForgetAfter(Instant const& t);
  void ForgetAfter(iterator it);

  // Deletes the trajectory points with a time in [begin, t[.  Preserves empty
  // segments and doesn't invalidate any segment iterator.
  void ForgetBefore(Instant const& t);
  void ForgetBefore(iterator it);

  // Return an error if downsampling was aborted.
  absl::Status Append(Instant const& t,
                      DegreesOfFreedom<Frame> const& degrees_of_freedom);

  // Merges |trajectory| (the source) into this object (the target).  The
  // operation processes pairs of segments taken from each trajectory and
  // proceeds as follows:
  // 1. If the source segment is empty (or missing), leave the target segment
  //    unchanged.
  // 2. If the target segment is empty (or missing), move the source segment
  //    into the target.  This includes its downsampling state.
  // 3. If both segments are nonempty, they must be nonoverlapping.  The points
  //    from the source segment are inserted in the target segment (possibly
  //    before the beginning of that segment).  The downsampling state of the
  //    result is that of the latest segment (the one with the largest times).
  void Merge(DiscreteTrajectory<Frame> trajectory);

  Instant t_min() const override;
  Instant t_max() const override;

  Position<Frame> EvaluatePosition(Instant const& t) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& t) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& t) const override;

  // The segments in |tracked| are restored at deserialization.  The points
  // denoted by |exact| are written and re-read exactly and are not affected by
  // any errors introduced by zfp compression.  The endpoints of each segment
  // are always exact.
  void WriteToMessage(
      not_null<serialization::DiscreteTrajectory*> message,
      std::vector<SegmentIterator> const& tracked,
      std::vector<iterator> const& exact) const;
  // Same as above, but only the points defined by [begin, end[ are written.
  void WriteToMessage(
      not_null<serialization::DiscreteTrajectory*> message,
      iterator begin,
      iterator end,
      std::vector<SegmentIterator> const& tracked,
      std::vector<iterator> const& exact) const;

  // |forks| must have a size appropriate for the |message| being deserialized
  // and the orders of the |forks| must be consistent during serialization and
  // deserialization.  All pointers designated by the pointers in |forks| must
  // be null at entry; they may be null at exit.
  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  static DiscreteTrajectory ReadFromMessage(
      serialization::DiscreteTrajectory const& message,
      std::vector<SegmentIterator*> const& tracked);

 private:
  using DownsamplingParameters =
      internal_discrete_trajectory_types::DownsamplingParameters;
  using Segments = internal_discrete_trajectory_types::Segments<Frame>;
  using SegmentByLeftEndpoint =
      absl::btree_map<Instant, typename Segments::iterator>;

  // This constructor leaves the list of segments empty (but allocated) as well
  // as the time-to-segment mapping.
  explicit DiscreteTrajectory(uninitialized_t);

  // Returns an iterator to a segment with extremities t1 and t2 such that
  // t ∈ [t1, t2[.  For the last segment, t2 is assumed to be +∞.  A 1-point
  // segment is never returned, unless it is the last one (because its upper
  // bound is assumed to be +∞).  Returns segment_by_left_endpoint_->end() iff
  // t is before the first time of the trajectory or if the trajectory is
  // empty().
  typename SegmentByLeftEndpoint::iterator FindSegment(Instant const& t);
  typename SegmentByLeftEndpoint::const_iterator
  FindSegment(Instant const& t) const;

  // Determines if this objects is in a consistent state, and returns an error
  // status with a relevant message if it isn't.
  absl::Status ConsistencyStatus() const;

  // Updates the segments self-pointers and the time-to-segment mapping after
  // segments have been spliced from |from| to |to|.  The iterator indicates the
  // segments to fix-up.
  static void AdjustAfterSplicing(
      DiscreteTrajectory& from,
      DiscreteTrajectory& to,
      typename Segments::iterator to_segments_begin);

  // Reads a pre-Hamilton downsampling message and return the downsampling
  // parameters and the start of the dense timeline.  The latter will have to be
  // converted to a number of points based on the deserialized timeline.
  static void ReadFromPreHamiltonMessage(
      serialization::DiscreteTrajectory::Downsampling const& message,
      DownsamplingParameters& downsampling_parameters,
      Instant& start_of_dense_timeline);

  // Reads a set of pre-Hamilton children.  Checks that there is only one child,
  // and that it is at the end of the preceding segment.  Append a segment to
  // the trajectory and returns an iterator to that segment.
  static SegmentIterator ReadFromPreHamiltonMessage(
      serialization::DiscreteTrajectory::Brood const& message,
      std::vector<SegmentIterator*> const& tracked,
      value_type const& fork_point,
      DiscreteTrajectory& trajectory);

  // Reads a pre-Hamilton trajectory, updating the tracked segments as needed.
  // If this is not the root of the trajectory, fork_point is set.
  static void ReadFromPreHamiltonMessage(
      serialization::DiscreteTrajectory const& message,
      std::vector<SegmentIterator*> const& tracked,
      std::optional<value_type> const& fork_point,
      DiscreteTrajectory& trajectory);

  // We need a level of indirection here to make sure that the pointer to
  // Segments in the DiscreteTrajectorySegmentIterator remain valid when the
  // DiscreteTrajectory moves.  This field is never null and never empty.
  not_null<std::unique_ptr<Segments>> segments_;

  // Maps time |t| to the last segment that start at time |t|.  Does not contain
  // entries for empty segments (at the beginning of the trajectory) or for
  // 1-point segments that are not the last at their time.  Empty iff the entire
  // trajectory is empty.  Always updated using |insert_or_assign| to override
  // any preexisting segment with the same endpoint.
  SegmentByLeftEndpoint segment_by_left_endpoint_;
};

}  // namespace internal_discrete_trajectory

using internal_discrete_trajectory::DiscreteTrajectory;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_body.hpp"
