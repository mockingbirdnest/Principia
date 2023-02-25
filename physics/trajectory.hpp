#pragma once

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {
namespace internal_trajectory {

using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using namespace principia::base::_not_null;

template<typename Frame>
class Trajectory {
 public:
  virtual ~Trajectory() = default;

  // The time range for which the trajectory can be evaluated is [t_min, t_max].
  // Note that it is not required that t_min ≤ t_max: for an empty trajectory,
  // t_min = +∞, and t_max = -∞.
  virtual Instant t_min() const = 0;
  virtual Instant t_max() const = 0;

  // Evaluates the trajectory at the given |time|, which must be in
  // [t_min(), t_max()].
  virtual Position<Frame> EvaluatePosition(Instant const& time) const = 0;
  virtual Velocity<Frame> EvaluateVelocity(Instant const& time) const = 0;
  virtual DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& time) const = 0;
};

}  // namespace internal_trajectory

using internal_trajectory::Trajectory;

}  // namespace physics
}  // namespace principia
