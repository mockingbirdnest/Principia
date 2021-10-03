#pragma once

#include <list>
#include <memory>

#include "absl/container/btree_map.h"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/trajectory.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {

template<typename Frame>
class DiscreteTrajectorySegment;

namespace internal_discrete_trajectory {

using base::not_null;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTrajectory2 : public Trajectory<Frame> {
 public:
  using Iterator = DiscreteTrajectoryIterator<Frame>;
  using SegmentIterator = DiscreteTrajectorySegmentIterator<Frame>;

  DiscreteTrajectory2() = default;

  // Moveable.
  DiscreteTrajectory2(DiscreteTrajectory2&&) = default;
  DiscreteTrajectory2& operator=(DiscreteTrajectory2&&) = default;
  DiscreteTrajectory2(const DiscreteTrajectory2&) = delete;
  DiscreteTrajectory2& operator=(const DiscreteTrajectory2&) = delete;

    Iterator begin() const;
  Iterator end() const;

  Iterator rbegin() const;
  Iterator rend() const;

  Iterator find(Instant const& t) const;

  Iterator lower_bound(Instant const& t) const;
  Iterator upper_bound(Instant const& t) const;

  SegmentIterator segments_begin() const;
  SegmentIterator segments_end() const;

  SegmentIterator segments_rbegin() const;
  SegmentIterator segments_rend() const;

  SegmentIterator NewSegment();

  DiscreteTrajectory DetachSegments(Iterator begin);
  SegmentIterator AttachSegments(
      DiscreteTrajectory&& trajectory);
  void DeleteSegments(Iterator begin);

  void ForgetAfter(Instant const& t);
  void ForgetAfter(Iterator begin);

  void ForgetBefore(Instant const& t);
  void ForgetBefore(Iterator end);

  void Append(Instant const& t,
              DegreesOfFreedom<Frame> const& degrees_of_freedom);

  Position<Frame> EvaluatePosition(Instant const& time) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& time) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& time) const override;

  void WriteToMessage(
      not_null<serialization::DiscreteTrajectory*> message,
      std::vector<SegmentIterator> const& tracked,
      std::vector<Iterator> const& exact) const;
  void WriteToMessage(
      not_null<serialization::DiscreteTrajectory*> message,
      Iterator begin, Iterator end,
      std::vector<SegmentIterator> const& tracked,
      std::vector<Iterator> const& exact) const;

  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  static not_null<std::unique_ptr<DiscreteTrajectory>> ReadFromMessage(
      serialization::DiscreteTrajectory const& message,
      std::vector<DiscreteTrajectory<Frame>**> const& tracked);

 private:
  using Segments = internal_discrete_trajectory_types::Segments<Frame>;

  Segments segments_;
};

}  // namespace internal_discrete_trajectory

template<typename Frame>
using DiscreteTrajectory2 = internal_discrete_trajectory::DiscreteTrajectory2;

}  // namespace principia
}  // namespace physics
