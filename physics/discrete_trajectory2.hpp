#pragma once

#include <iterator>
#include <list>
#include <memory>
#include <vector>

#include "absl/container/btree_map.h"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
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

  iterator find(Instant const& t) const;

  iterator lower_bound(Instant const& t) const;
  iterator upper_bound(Instant const& t) const;

  SegmentRange segments() const;
  // TODO(phl): In C++20 this should be a reverse_view on segments.
  ReverseSegmentRange rsegments() const;

  SegmentIterator NewSegment();

  DiscreteTrajectory2 DetachSegments(iterator begin);
  SegmentIterator AttachSegments(DiscreteTrajectory2&& trajectory);
  void DeleteSegments(iterator begin);

  void ForgetAfter(Instant const& t);
  void ForgetAfter(iterator begin);

  void ForgetBefore(Instant const& t);
  void ForgetBefore(iterator end);

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
  void WriteToMessage(
      not_null<serialization::DiscreteTrajectory*> message,
      iterator begin, iterator end,
      std::vector<SegmentIterator> const& tracked,
      std::vector<iterator> const& exact) const;

  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  static not_null<std::unique_ptr<DiscreteTrajectory2>> ReadFromMessage(
      serialization::DiscreteTrajectory const& message,
      std::vector<DiscreteTrajectory2**> const& tracked);

 private:
  using Segments = internal_discrete_trajectory_types::Segments<Frame>;

  DiscreteTrajectorySegment<Frame>& FindSegment(Instant const& t);
  DiscreteTrajectorySegment<Frame> const& FindSegment(Instant const& t) const;

  // We need a level of indirection here to make sure that the pointer to
  // Segments in the DiscreteTrajectorySegmentIterator remain valid when the
  // DiscreteTrajectory moves.
  //TODO(phl):never empty
  not_null<std::unique_ptr<Segments>> segments_;

  //TODO(phl): Use --upper_bound(t) to access, check for t_min/t_max;
  absl::btree_map<Instant, Segments::const_iterator> segment_by_left_endpoint_;
};

}  // namespace internal_discrete_trajectory2

using internal_discrete_trajectory2::DiscreteTrajectory2;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory2_body.hpp"
