#pragma once

#include <iterator>
#include <list>
#include <memory>
#include <ranges>
#include <vector>

#include "absl/container/btree_map.h"
#include "absl/status/status.h"
#include "base/concepts.hpp"
#include "base/not_null.hpp"
#include "base/tags.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/trajectory.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace _discrete_trajectory {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_not_null;
using namespace principia::base::_tags;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory_iterator;
using namespace principia::physics::_discrete_trajectory_segment;
using namespace principia::physics::_discrete_trajectory_segment_iterator;
using namespace principia::physics::_discrete_trajectory_types;
using namespace principia::physics::_trajectory;

//TODO(phl)comments
template<typename Frame>
class DiscreteTrajectoryView : public Trajectory<Frame> {
 public:
  using key_type = typename Timeline<Frame>::key_type;
  using value_type = typename Timeline<Frame>::value_type;

  using iterator = DiscreteTrajectoryIterator<Frame>;
  using const_iterator = DiscreteTrajectoryIterator<Frame>;
  using reference = value_type const&;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using SegmentIterator = DiscreteTrajectorySegmentIterator<Frame>;
  using ReverseSegmentIterator = std::reverse_iterator<SegmentIterator>;
  using SegmentRange = std::ranges::subrange<SegmentIterator,
                                             SegmentIterator,
                                             std::ranges::subrange_kind::sized>;

  DiscreteTrajectoryView(DiscreteTrajectory<Frame> const& trajectory,
                         const_iterator begin,
                         const_iterator end);
  //TODO(phl)constructors

  // Copyable.
  DiscreteTrajectoryView(DiscreteTrajectoryView&&) = default;
  DiscreteTrajectoryView& operator=(DiscreteTrajectoryView&&) = default;
  DiscreteTrajectoryView(const DiscreteTrajectoryView&) = delete;
  DiscreteTrajectoryView& operator=(const DiscreteTrajectoryView&) = delete;

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
  std::ranges::reverse_view<SegmentRange> rsegments() const;

  Instant t_min() const override;
  Instant t_max() const override;

  Position<Frame> EvaluatePosition(Instant const& t) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& t) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& t) const override;

  // Not serializable.

 private:
  not_null<DiscreteTrajectory<Frame>*> trajectory_;
  const_iterator begin_;
  const_iterator end_;
  Instant t_min;
  Instant t_max;
};

}  // namespace internal

using internal::DiscreteTrajectory;

}  // namespace _discrete_trajectory
}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_view_body.hpp"
