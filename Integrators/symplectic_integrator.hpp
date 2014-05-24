#pragma once

#include <vector>

namespace principia {
namespace integrators {

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

  // TODO(phl): See if making these template parameters would improve
  // performance.
  typedef void (*AutonomousRightHandSideComputation)(
                     std::vector<double> const& y,
                     std::vector<double>* result);
  typedef void (*RightHandSideComputation)(double const t,
                                           std::vector<double> const& y,
                                           std::vector<double>* result);

  // TODO(phl): Rework the struct names, maybe promote them to classes.
  struct Parameters {
    Parameters();
    // The initial positions for each dimension.
    std::vector<double> q0;
    // The initial momenta for each dimension.
    std::vector<double> p0;
    // The initial position errors for each dimension, null if no error.
    std::vector<double>* q_error;
    // The initial momentum errors for each dimension, null if no error.
    std::vector<double>* p_error;
    // The starting time of the resolution.
    double t0;
    // The ending time of the resolution.
    double tmax;
    // The time step.
    double Δt;
    // The error on ????
    double t_error;
    // The coefficients of the integrator.
    std::vector<std::vector<double>> coefficients;
    // The sampling ????
    int sampling_period;
  };

  struct TimeseriesAndError {
    // Indexed by time step.
    std::vector<double> quantities;
    double error;
  };

  struct Solution {
    // Indexed by dimension.
    std::vector<TimeseriesAndError> momentum;
    std::vector<TimeseriesAndError> position;
    TimeseriesAndError time;
  };

  // Takes ownership of the pointers in |parameters|.
  // TODO(phl): Pass the function pointers at construction?
  virtual void Solve(RightHandSideComputation const compute_force,
                     AutonomousRightHandSideComputation const compute_velocity,
                     Parameters const& parameters,
                     Solution* solution) = 0;
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_integrator_body.hpp"
