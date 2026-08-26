#pragma once

#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/trajectory.hpp"

namespace principia {
namespace physics {
namespace _trajectory_view {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_trajectory;

// A read-only view of a range of a `Trajectory`.  This class is copyable.
template<typename Frame>
class TrajectoryView : public Trajectory<Frame> {
 public:
  // The `Trajectory` passed at construction must outlive this object.

  // A convenience constructor for building a view that covers an entire
  // trajectory.
  explicit TrajectoryView(not_null<Trajectory<Frame> const*> trajectory);

  // Constructs a view that covers the parts of the `trajectory` whose time is
  // in the range [`t_min`, `t_max`].  If `t_max < t_min`, the view is empty and
  // its `t_min()` and `t_max()` are +∞ and -∞, respectively.  Otherwise `t_min`
  // and `t_max` must be in the time range of the `trajectory`.
  TrajectoryView(not_null<Trajectory<Frame> const*> trajectory,
                 Instant const& t_min,
                 Instant const& t_max);

  bool empty() const;

  // These functions return +∞ and -∞, respectively, for an empty view.
  Instant t_min() const override;
  Instant t_max() const override;

  // Modifies this view by intersecting its time range with the given
  // [`t_min`, `t_max`].  This may result in an empty view.
  virtual void Restrict(Instant const& t_min, Instant const& t_max);

  // `t` must be in the range [`this->t_min()`, `this->t_max()`].
  Position<Frame> EvaluatePosition(Instant const& t) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& t) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& t) const override;

  // Not serializable.

 protected:
  Trajectory<Frame> const& trajectory() const;

 private:
  not_null<Trajectory<Frame> const*> trajectory_;
  Instant t_min_;
  Instant t_max_;
};

template<typename Frame>
TrajectoryView(Trajectory<Frame> const*) -> TrajectoryView<Frame>;

template<typename Frame>
TrajectoryView(Trajectory<Frame> const*, Instant const&, Instant const&)
    -> TrajectoryView<Frame>;

}  // namespace internal

using internal::TrajectoryView;

}  // namespace _trajectory_view
}  // namespace physics
}  // namespace principia

#include "physics/trajectory_view_body.hpp"
