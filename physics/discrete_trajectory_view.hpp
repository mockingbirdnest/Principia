#pragma once

#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/trajectory.hpp"

namespace principia {
namespace physics {
namespace _discrete_trajectory_view {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_discrete_trajectory_iterator;
using namespace principia::physics::_discrete_trajectory_types;
using namespace principia::physics::_trajectory;

// A view of a range of a `DiscreteTrajectory`.  This class is copyable.
template<typename Frame>
class DiscreteTrajectoryView : public Trajectory<Frame> {
 public:
  using key_type = typename Timeline<Frame>::key_type;
  using value_type = typename Timeline<Frame>::value_type;

  using iterator = DiscreteTrajectoryIterator<Frame>;
  using const_iterator = DiscreteTrajectoryIterator<Frame>;
  using reference = value_type const&;
  using reverse_iterator = std::reverse_iterator<iterator>;

  // The `DiscreteTrajectory` passed at construction must outlive this object.

  // A convenience constructor for building a view that covers an entire
  // trajectory.
  explicit DiscreteTrajectoryView(
      not_null<DiscreteTrajectory<Frame> const*> trajectory);

  // Constructs a view that covers the range [`begin`, `end`[.
  DiscreteTrajectoryView(not_null<DiscreteTrajectory<Frame> const*> trajectory,
                         const_iterator begin,
                         const_iterator end);

  // Constructs a view that covers all the points in the `trajectory` whose time
  // is in the range [`t_min`, `t_max`].  Note that the resulting view `v` has
  // `t_min <= v.t_min()` and `v.t_max() <= t_max`.  It may be empty if there
  // is not point in the given time interval.
  DiscreteTrajectoryView(not_null<DiscreteTrajectory<Frame> const*> trajectory,
                         Instant const& t_min,
                         Instant const& t_max);

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

  // No `segments` or `rsegments` as that would expose parts of the trajectory
  // outside of the range of the view.

  // These functions return +∞ and -∞, respectively, for an empty view.
  Instant t_min() const override;
  Instant t_max() const override;

  // `t` must be in the range [`this->t_min()`, `this->t_max()`].
  Position<Frame> EvaluatePosition(Instant const& t) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& t) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& t) const override;

  // Not serializable.

 private:
  not_null<DiscreteTrajectory<Frame> const*> trajectory_;
  const_iterator begin_;
  const_iterator end_;
  Instant t_min_;
  Instant t_max_;
};

}  // namespace internal

using internal::DiscreteTrajectoryView;

}  // namespace _discrete_trajectory_view
}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_view_body.hpp"
