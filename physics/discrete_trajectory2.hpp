#pragma once

#include <iterator>
#include <list>
#include <memory>
#include <vector>

#include "absl/container/btree_map.h"
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

namespace internal_discrete_trajectory2 {

using base::not_null;
using base::uninitialized_t;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTrajectory2 : public Trajectory<Frame> {
 public:
  using key_type =
      typename internal_discrete_trajectory_types::Timeline<Frame>::key_type;
  using value_type =
      typename internal_discrete_trajectory_types::Timeline<Frame>::value_type;

  using iterator = DiscreteTrajectoryIterator<Frame>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using SegmentIterator = DiscreteTrajectorySegmentIterator<Frame>;
  using ReverseSegmentIterator = std::reverse_iterator<SegmentIterator>;
  using SegmentRange = DiscreteTrajectorySegmentRange<SegmentIterator>;
  using ReverseSegmentRange =
      DiscreteTrajectorySegmentRange<ReverseSegmentIterator>;

  DiscreteTrajectory2();

  // Moveable.
  DiscreteTrajectory2(DiscreteTrajectory2&&) = default;
  DiscreteTrajectory2& operator=(DiscreteTrajectory2&&) = default;
  DiscreteTrajectory2(const DiscreteTrajectory2&) = delete;
  DiscreteTrajectory2& operator=(const DiscreteTrajectory2&) = delete;

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

  DiscreteTrajectory2 DetachSegments(SegmentIterator begin);
  SegmentIterator AttachSegments(DiscreteTrajectory2&& trajectory);
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
  static DiscreteTrajectory2 ReadFromMessage(
      serialization::DiscreteTrajectory const& message,
      std::vector<SegmentIterator*> const& tracked);

 private:
  using DownsamplingParameters =
      internal_discrete_trajectory_types::DownsamplingParameters;
  using Segments = internal_discrete_trajectory_types::Segments<Frame>;

  // This constructor leaves the list of segments empty (but allocated) as well
  // as the time-to-segment mapping.
  explicit DiscreteTrajectory2(uninitialized_t);

  typename Segments::iterator FindSegment(Instant const& t);
  typename Segments::const_iterator FindSegment(Instant const& t) const;

  // Updates the segments self-pointers and the time-to-segment mapping after
  // segments have been spliced from |from| to |to|.  The iterators indicate the
  // segments to fix-up.
  static void AdjustAfterSplicing(
      DiscreteTrajectory2& from,
      DiscreteTrajectory2& to,
      typename Segments::iterator to_segments_begin,
      std::reverse_iterator<typename Segments::iterator> to_segments_rend);

  //TODO(phl):comment
  static void ReadFromPreΖήνωνMessage(
      serialization::DiscreteTrajectory::Downsampling const& message,
      DownsamplingParameters& downsampling_parameters,
      Instant& start_of_dense_timeline);

  static SegmentIterator ReadFromPreΖήνωνMessage(
      serialization::DiscreteTrajectory::Brood const& message,
      std::vector<SegmentIterator*> const& tracked,
      value_type const& fork_point,
      DiscreteTrajectory2& trajectory);

  static void ReadFromPreΖήνωνMessage(
      serialization::DiscreteTrajectory const& message,
      std::vector<SegmentIterator*> const& tracked,
      DiscreteTrajectory2& trajectory);

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

}  // namespace internal_discrete_trajectory2

using internal_discrete_trajectory2::DiscreteTrajectory2;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory2_body.hpp"
