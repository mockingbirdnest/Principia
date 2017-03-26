
#pragma once

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {
namespace internal_trajectory {

using base::not_null;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;

template<typename Frame>
class Trajectory {
 public:
  virtual ~Trajectory() = default;

  // Derived classes of |Trajectory| should also derive |Hint| as appropriate.
  class Hint {
   public:
    virtual ~Hint() = default;
   protected:
    Hint() = default;
  };

  // The time range for which the trajectory can be evaluated.  Note that it is
  // not required that t_min ≤ t_max: for an empty trajectory, t_min = +∞, and
  // t_max = -∞.
  virtual Instant t_min() const = 0;
  virtual Instant t_max() const = 0;

  // Gets a |Hint| of the appropriate derived type.
  virtual not_null<std::unique_ptr<Hint>> GetHint() const = 0;

  // Evaluates the trajectory at the given |time|, which must be in
  // [t_min(), t_max()].  The |hint| may be used to speed up evaluation in
  // increasing time order.  It may be a nullptr (in which case no speed-up
  // takes place), otherwise it must be of the derived |Hint| type corresponding
  // to the |Trajectory|.
  virtual Position<Frame> EvaluatePosition(Instant const& time,
                                           Hint* hint) const = 0;
  virtual Velocity<Frame> EvaluateVelocity(Instant const& time,
                                           Hint* hint) const = 0;
  virtual DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& time,
      Hint* hint) const = 0;
};

}  // namespace internal_trajectory

using internal_trajectory::Trajectory;

}  // namespace physics
}  // namespace principia
