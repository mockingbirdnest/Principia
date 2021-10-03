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
  DiscreteTrajectory2() = default;

  DiscreteTrajectoryIterator begin() const;
  DiscreteTrajectoryIterator end() const;

  DiscreteTrajectoryIterator rbegin() const;
  DiscreteTrajectoryIterator rend() const;

  DiscreteTrajectoryIterator find(Instant const& t) const;

  DiscreteTrajectoryIterator lower_bound(Instant const& t) const;
  DiscreteTrajectoryIterator upper_bound(Instant const& t) const;

  DiscreteTrajectorySegmentIterator segments_begin() const;
  DiscreteTrajectorySegmentIterator segments_end() const;

  DiscreteTrajectorySegmentIterator segments_rbegin() const;
  DiscreteTrajectorySegmentIterator segments_rend() const;

  DiscreteTrajectorySegmentIterator NewSegment();

  DiscreteTrajectory DetachSegments(DiscreteTrajectoryIterator begin);
  DiscreteTrajectorySegmentIterator AttachSegments(
      DiscreteTrajectory&& trajectory);
  void DeleteSegments(DiscreteTrajectoryIterator begin);

  void ForgetAfter(Instant const& t);
  void ForgetAfter(DiscreteTrajectoryIterator begin);

  void ForgetBefore(Instant const& t);
  void ForgetBefore(DiscreteTrajectoryIterator end);

  void Append(Instant const& t,
              DegreesOfFreedom<Frame> const& degrees_of_freedom);

  Position<Frame> EvaluatePosition(Instant const& time) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& time) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& time) const override;

  void WriteToMessage(
      not_null<serialization::DiscreteTrajectory*> message,
      std::vector<DiscreteTrajectorySegmentIterator> const& tracked,
      std::vector<DiscreteTrajectoryIterator> const& exact) const;
  void WriteToMessage(
      not_null<serialization::DiscreteTrajectory*> message,
      DiscreteTrajectoryIterator begin,
      DiscreteTrajectoryIterator end,
      std::vector<DiscreteTrajectorySegmentIterator> const& tracked,
      std::vector<DiscreteTrajectoryIterator> const& exact) const;

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
