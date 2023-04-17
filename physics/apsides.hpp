#pragma once

#include <functional>

#include "absl/status/status.h"
#include "base/constant_function.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/trajectory.hpp"

namespace principia {
namespace physics {
namespace _apsides {
namespace internal {

using namespace principia::base::_constant_function;
using namespace principia::geometry::_grassmann;

// Computes the apsides with respect to |reference| for the section given by
// |begin| and |end| of |trajectory|.  Appends to the given output trajectories
// one point for each apsis.
template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& reference,
                    Trajectory<Frame> const& trajectory,
                    typename DiscreteTrajectory<Frame>::iterator begin,
                    typename DiscreteTrajectory<Frame>::iterator end,
                    int max_points,
                    DiscreteTrajectory<Frame>& apoapsides,
                    DiscreteTrajectory<Frame>& periapsides);

// Computes the crossings of the section given by |begin| and |end| of
// |trajectory| with the xy plane.  Appends the crossings that go towards the
// |north| side of the xy plane to |ascending|, and those that go away from the
// |north| side to |descending|.
// Nodes for which |predicate| returns false are excluded.
template<typename Frame, typename Predicate = ConstantFunction<bool>>
absl::Status ComputeNodes(Trajectory<Frame> const& trajectory,
                          typename DiscreteTrajectory<Frame>::iterator begin,
                          typename DiscreteTrajectory<Frame>::iterator end,
                          Vector<double, Frame> const& north,
                          int max_points,
                          DiscreteTrajectory<Frame>& ascending,
                          DiscreteTrajectory<Frame>& descending,
                          Predicate predicate = Identically(true));

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

}  // namespace internal

using internal::ComputeApsides;
using internal::ComputeNodes;

}  // namespace _apsides
}  // namespace physics
}  // namespace principia

#include "physics/apsides_body.hpp"
