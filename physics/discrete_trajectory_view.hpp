#pragma once

#include <ranges>

#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/trajectory_view.hpp"

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
using namespace principia::physics::_trajectory_view;

// A read-only view of a range of a `DiscreteTrajectory`.  This class is
// copyable.
template<typename Frame>
class DiscreteTrajectoryView : public TrajectoryView<Frame> {
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
  // is in the interval [`t_min`, `t_max`].  Note that iterating over the
  // (discrete) points of the trajectory may start with a point above `t_min`
  // and end with a point below `t_max`.  The view may even be empty if there is
  // no point in the given time interval.  But it is still possible to evaluate
  // the degrees of freedom over the entire interval [`t_min`, `t_max`].
  DiscreteTrajectoryView(not_null<DiscreteTrajectory<Frame> const*> trajectory,
                         Instant const& t_min,
                         Instant const& t_max);

  reference front() const;
  reference back() const;

  iterator begin() const;
  iterator end() const;

  reverse_iterator rbegin() const;
  reverse_iterator rend() const;

  // Beware: `empty()` implies `size() == 0` but not the reverse.  `empty()`
  // means that there exist no valid times to evaluate the degrees of freedom of
  // the view.  `size() == 0` means that iterating over the (discrete) points of
  // the view yields no point.  A trajectory constructed with a time interval
  // that is between two points of the discrete trajectory is not empty but its
  // size is 0.
  std::int64_t size() const;

  iterator find(Instant const& t) const;

  iterator lower_bound(Instant const& t) const;
  iterator upper_bound(Instant const& t) const;

  // No `segments` or `rsegments` as that would expose parts of the trajectory
  // outside of the range of the view.

  // Modifies this view by intersecting its time range with the given
  // [`t_min`, `t_max`].  This may result in an empty view or a view of size 0.
  void Restrict(Instant const& t_min, Instant const& t_max) override;

  // Same as above, but with iterators.
  void Restrict(const_iterator begin, const_iterator end);

  // Not serializable.

 private:
  DiscreteTrajectory<Frame> const& trajectory() const;

  static Instant TrajectoryViewTMin(DiscreteTrajectory<Frame> const& trajectory,
                                    const_iterator const begin,
                                    const_iterator const end);
  static Instant TrajectoryViewTMax(DiscreteTrajectory<Frame> const& trajectory,
                                    const_iterator const begin,
                                    const_iterator const end);

  const_iterator begin_;
  const_iterator end_;
};

template<typename Frame>
DiscreteTrajectoryView(DiscreteTrajectory<Frame> const*)
    -> DiscreteTrajectoryView<Frame>;

template<typename Frame>
DiscreteTrajectoryView(DiscreteTrajectory<Frame> const*,
                       typename DiscreteTrajectoryView<Frame>::const_iterator,
                       typename DiscreteTrajectoryView<Frame>::const_iterator)
    -> DiscreteTrajectoryView<Frame>;

template<typename Frame>
DiscreteTrajectoryView(DiscreteTrajectory<Frame> const*,
                       Instant const&,
                       Instant const&)
    -> DiscreteTrajectoryView<Frame>;

}  // namespace internal

using internal::DiscreteTrajectoryView;

}  // namespace _discrete_trajectory_view
}  // namespace physics
}  // namespace principia

namespace std {
namespace ranges {
template<typename Frame>
inline constexpr bool enable_borrowed_range<
    principia::physics::_discrete_trajectory_view::DiscreteTrajectoryView<
        Frame>> = true;
}  // namespace ranges
}  // namespace std

#include "physics/discrete_trajectory_view_body.hpp"
