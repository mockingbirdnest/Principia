#pragma once

#include <functional>

#include "physics/discrete_trajectory.hpp"
#include "physics/trajectory.hpp"

namespace principia {
namespace physics {
namespace internal_apsides {

using geometry::Vector;

// Computes the apsides with respect to |reference| for the discrete trajectory
// segment given by |begin| and |end|.  Appends to the given trajectories one
// point for each apsis.
template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& reference,
                    typename DiscreteTrajectory<Frame>::Iterator begin,
                    typename DiscreteTrajectory<Frame>::Iterator end,
                    DiscreteTrajectory<Frame>& apoapsides,
                    DiscreteTrajectory<Frame>& periapsides);

// REMOVE BEFORE FLIGHT: move to base.
template<typename T>
struct constant_function {
  template<typename... Args>
  constexpr T operator()(Args&&...) {
    return value;
  }

  T value;
};

template <typename T>
constexpr constant_function<T> Everywhere(T value) {
  return {value};
}

// Computes the crossings of the discrete trajectory segment given by |begin|
// and |end| with the xy plane.  Appends the crossings that go towards the
// |north| side of the xy plane to |ascending|, and those that go away from the
// |north| side to |descending|.
/// Nodes for which |filter| returns false are excluded.
template<typename Frame, typename Predicate = constant_function<bool>>
void ComputeNodes(typename DiscreteTrajectory<Frame>::Iterator begin,
                  typename DiscreteTrajectory<Frame>::Iterator end,
                  Vector<double, Frame> const& north,
                  DiscreteTrajectory<Frame>& ascending,
                  DiscreteTrajectory<Frame>& descending,
                  Predicate filter = Everywhere(true));

// TODO(egg): when we can usefully iterate over an arbitrary |Trajectory|, move
// the following from |Ephemeris|.
#if 0
template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& trajectory1,
                    Trajectory<Frame> const& trajectory2,
                    DiscreteTrajectory<Frame>& apoapsides1,
                    DiscreteTrajectory<Frame>& periapsides1,
                    DiscreteTrajectory<Frame>& apoapsides2,
                    DiscreteTrajectory<Frame>& periapsides2);
#endif

}  // namespace internal_apsides

using internal_apsides::ComputeApsides;
using internal_apsides::ComputeNodes;

}  // namespace physics
}  // namespace principia

#include "physics/apsides_body.hpp"
