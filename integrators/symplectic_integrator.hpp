#pragma once

#include <vector>

#include "quantities/quantities.hpp"

using principia::quantities::Time;

namespace principia {
namespace integrators {

// A simple container for a scalar value and the related error.  The
// constructor is not explicit to make it easy to construct an object with no
// error.
template<typename Scalar>
struct DoublePrecision {
  DoublePrecision() = default;
  DoublePrecision(Scalar const& value);  // NOLINT(runtime/explicit)
  Scalar value;
  Scalar error;
};

template<typename Position, typename Momentum>
class SymplecticIntegrator {
 public:
  SymplecticIntegrator() = default;
  virtual ~SymplecticIntegrator() = default;

  // The coefficients of the integrator.
  using Coefficients = std::vector<std::vector<double>>;

  // The entire state of the system at a given time.  The vectors are indexed by
  // dimension.
  struct SystemState {
    std::vector<DoublePrecision<Position>> positions;
    std::vector<DoublePrecision<Momentum>> momenta;
    DoublePrecision<Time> time;
  };

  // TODO(phl): Should we remove this class entirely and just pass 4 parameters
  // to Solve?
  struct Parameters {
    // The initial state of the system.
    SystemState initial;
    // The ending time of the resolution.
    Time tmax;
    // The time step.
    Time Δt;
    // To save memory, we only return a datapoint every sampling_period steps
    // (for trajectory plotting), as well as the result from the last step. If
    // sampling_period == 0, we only return the result from the last step
    // (that's for when we just want to advance the system, not to plot its
    // evolution).
    // NOTE(eggrobin): The images in the OP of the forum thread show the problem
    // with the current approach: with reasonable sampling_periods, the plotted
    // trajectory sometimes becomes polygonal at high velocities, while points
    // are wasted at low velocities. At some point I think this should be
    // handled with a function that evaluates the velocity in the plot frame to
    // decide when to sample.  Plotting some sort of higher-order spline, rather
    // than a polygon, would help, but isn't enough.
    int sampling_period;
  };

  // Initialize the integrator with the given |coefficients|.  Must be called
  // before calling Solve.
  virtual void Initialize(Coefficients const& coefficients) = 0;

  // NOTE(phl): This is part of the contract of SymplecticIntegrator but it
  // cannot be written in C++ because a template cannot be virtual.
  //
  // Takes ownership of the pointers in |parameters|.
  // template<typename AutonomousRightHandSideComputation,
  //          typename RightHandSideComputation>
  // void Solve(RightHandSideComputation const compute_force,
  //            AutonomousRightHandSideComputation const compute_velocity,
  //            Parameters const& parameters,
  //            std::vector<SystemState>* solution);
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_integrator_body.hpp"
