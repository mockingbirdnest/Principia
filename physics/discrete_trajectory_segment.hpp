#pragma once

#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment {

using geometry::Instant;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTrajectorySegment {
 public:
  DiscreteTrajectorySegment() = default;

  DiscreteTrajectoryIterator<Frame> begin() const;
  DiscreteTrajectoryIterator<Frame> end() const;

  DiscreteTrajectoryIterator<Frame> rbegin() const;
  DiscreteTrajectoryIterator<Frame> rend() const;

  DiscreteTrajectoryIterator<Frame> find(Instant const& t) const;

  DiscreteTrajectoryIterator<Frame> lower_bound(Instant const& t) const;
  DiscreteTrajectoryIterator<Frame> upper_bound(Instant const& t) const;

 private:
  using Timeline = internal_discrete_trajectory_types::Timeline<Frame>;

  void Append(Instant const& t,
              DegreesOfFreedom<Frame> const& degrees_of_freedom);

  void ForgetAfter(Instant const& t);
  void ForgetAfter(Timeline::const_iterator begin);

  void ForgetBefore(Instant const& t);
  void ForgetBefore(Timeline::const_iterator end);

  DiscreteTrajectorySegmentIterator<Frame> that_;

  Timeline timeline_;
  absl::btree_set<Instant> dense_points_;
};

}  // namespace internal_discrete_trajectory_segment

template<typename Frame>
using DiscreteTrajectorySegment =
    internal_discrete_trajectory_segment::DiscreteTrajectorySegment;

}  // namespace principia
}  // namespace physics
