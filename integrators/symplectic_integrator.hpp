#pragma once

#include <vector>

#include "quantities/quantities.hpp"

using principia::quantities::Time;

namespace principia {
namespace integrators {

template<typename Position, typename Momentum>
class SymplecticIntegrator {
 public:
  SymplecticIntegrator() = default;
  virtual ~SymplecticIntegrator() = default;

  // TODO(phl): I'd prefer this API to be typed with the untyping done
  // internally.
  // NOTE(eggrobin): This namespace is supposed to contain general-purpose
  // integrators. If you were to type it, you'd have to define
  // Product<double, double>, Product<float, float> etc. in the same place as
  // the other products (quantities would hardly be appropriate in that case).

  // The coefficients of the integrator.
  typedef std::vector<std::vector<double>> Coefficients;

  // TODO(phl): Rework the struct names, maybe promote them to classes.
  struct Parameters {
    Parameters();
    // The initial positions for each dimension.
    std::vector<Position> q0;
    // The initial momenta for each dimension.
    std::vector<Momentum> p0;
    // The initial position errors for each dimension, null if no error.
    std::vector<Position>* q_error;
    // The initial momentum errors for each dimension, null if no error.
    std::vector<Momentum>* p_error;
    // The starting time of the resolution.
    Time t0;
    // The ending time of the resolution.
    Time tmax;
    // The time step.
    Time Δt;
    // The error on the starting time.
    Time t_error;
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

  template<typename Scalar>
  struct TimeseriesAndError {
    // Indexed by time step.
    std::vector<Scalar> quantities;
    Scalar error;
  };

  struct Solution {
    // Indexed by dimension.
    std::vector<TimeseriesAndError<Position>> position;
    std::vector<TimeseriesAndError<Momentum>> momentum;
    TimeseriesAndError<Time> time;
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
  //            Solution* solution);
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_integrator_body.hpp"
