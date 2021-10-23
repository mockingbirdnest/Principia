﻿#pragma once

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

namespace internal_discrete_traject0ry {

using base::not_null;
using base::uninitialized_t;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTraject0ry : public Trajectory<Frame> {
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

  DiscreteTraject0ry();

  // Moveable.
  DiscreteTraject0ry(DiscreteTraject0ry&&) = default;
  DiscreteTraject0ry& operator=(DiscreteTraject0ry&&) = default;
  DiscreteTraject0ry(const DiscreteTraject0ry&) = delete;
  DiscreteTraject0ry& operator=(const DiscreteTraject0ry&) = delete;

  reference front() const;
  reference back() const;

  iterator begin() const;
  iterator end() const;

  reverse_iterator rbegin() const;
  reverse_iterator rend() const;

  bool empty() const;
  std::int64_t size() const;

  iterator find(Instant const& t) const;

  iterator lower_bound(Instant const& t) const;
  iterator upper_bound(Instant const& t) const;

  SegmentRange segments() const;
  // TODO(phl): In C++20 this should be a reverse_view on segments.
  ReverseSegmentRange rsegments() const;

  SegmentIterator NewSegment();

  DiscreteTraject0ry DetachSegments(SegmentIterator begin);
  SegmentIterator AttachSegments(DiscreteTraject0ry&& trajectory);
  void DeleteSegments(SegmentIterator begin);

  void ForgetAfter(Instant const& t);
  void ForgetBefore(Instant const& t);

  void Append(Instant const& t,
              DegreesOfFreedom<Frame> const& degrees_of_freedom);

  Instant t_min() const override;
  Instant t_max() const override;

  Position<Frame> EvaluatePosition(Instant const& t) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& t) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::DiscreteTrajectory*> message,
      std::vector<SegmentIterator> const& tracked,
      std::vector<iterator> const& exact) const;
  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  static DiscreteTraject0ry ReadFromMessage(
      serialization::DiscreteTrajectory const& message,
      std::vector<SegmentIterator*> const& tracked);

 private:
  using DownsamplingParameters =
      internal_discrete_trajectory_types::DownsamplingParameters;
  using Segments = internal_discrete_trajectory_types::Segments<Frame>;

  // This constructor leaves the list of segments empty (but allocated) as well
  // as the time-to-segment mapping.
  explicit DiscreteTraject0ry(uninitialized_t);

  typename Segments::iterator FindSegment(Instant const& t);
  typename Segments::const_iterator FindSegment(Instant const& t) const;

  // Checks if this objects is in a consistent state, and returns an error
  // status with a relevant message if it isn't.
  absl::Status ValidateConsistency() const;

  // Updates the segments self-pointers and the time-to-segment mapping after
  // segments have been spliced from |from| to |to|.  The iterators indicate the
  // segments to fix-up.
  static void AdjustAfterSplicing(
      DiscreteTraject0ry& from,
      DiscreteTraject0ry& to,
      typename Segments::iterator to_segments_begin,
      std::reverse_iterator<typename Segments::iterator> to_segments_rend);

  // Reads a pre-Ζήνων downsampling message and return the downsampling
  // parameters and the start of the dense timeline.  The latter will have to be
  // converted to a number of points based on the deserialized timeline.
  static void ReadFromPreΖήνωνMessage(
      serialization::DiscreteTrajectory::Downsampling const& message,
      DownsamplingParameters& downsampling_parameters,
      Instant& start_of_dense_timeline);

  // Reads a set of pre-Ζήνων children.  Checks that there is only one child,
  // and that it is at the end of the preceding segment.  Append a segment to
  // the trajectory and returns an iterator to that segment.
  static SegmentIterator ReadFromPreΖήνωνMessage(
      serialization::DiscreteTrajectory::Brood const& message,
      std::vector<SegmentIterator*> const& tracked,
      value_type const& fork_point,
      DiscreteTraject0ry& trajectory);

  // Reads a pre-Ζήνων trajectory, updating the tracked segments as needed.  If
  // this is not the root of the trajectory, fork_point is set.
  static void ReadFromPreΖήνωνMessage(
      serialization::DiscreteTrajectory const& message,
      std::vector<SegmentIterator*> const& tracked,
      std::optional<value_type> const& fork_point,
      DiscreteTraject0ry& trajectory);

  // We need a level of indirection here to make sure that the pointer to
  // Segments in the DiscreteTrajectorySegmentIterator remain valid when the
  // DiscreteTrajectory moves.  This field is never null and never empty.
  not_null<std::unique_ptr<Segments>> segments_;

  // This list is never empty.  For an empty trajectory, there is a sentinel
  // with time -∞ denoting the single segment of the trajectory.  As soon as a
  // point is appended to the trajectory, the sentinel is removed and a bona
  // fide entry replaces it.  To access the segment for time t, use
  // |--upper_bound(t)|.
  absl::btree_map<Instant,
                  typename Segments::iterator> segment_by_left_endpoint_;
};

}  // namespace internal_discrete_traject0ry

using internal_discrete_traject0ry::DiscreteTraject0ry;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_traject0ry_body.hpp"