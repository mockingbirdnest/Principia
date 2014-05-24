#pragma once

#include <vector>

namespace principia {
namespace integrators {

class Integrator {
 public:
  Integrator();
  virtual ~Integrator();

  // TODO(phl): I'd prefer this API to be typed with the untyping done
  // internally.

  typedef void (*AutonomousRightHandSideComputation)(
                     std::vector<double> const& y,
                     std::vector<double>* result);
  typedef void (*RightHandSideComputation)(double const t,
                                           std::vector<double> const& y,
                                           std::vector<double>* result);

  // TODO(phl): Rework the struct names, maybe promote them to classes.
  struct Parameters {
    Parameters();
    std::vector<double> q0;
    std::vector<double> p0;
    std::vector<double>* q_error;
    std::vector<double>* p_error;
    double t0;
    double tmax;
    double Δt;
    double t_error;
    std::vector<std::vector<double>> coefficients;
    int sampling_period;
  };

  struct QuantitiesAndError {
    std::vector<double> quantities;
    double error;
  };

  struct Solution {
    std::vector<QuantitiesAndError> momentum;
    std::vector<QuantitiesAndError> position;
    QuantitiesAndError time;
  };

  // Takes ownership of the pointers in |parameters|.
  virtual void Increment(
      RightHandSideComputation const compute_force,
      AutonomousRightHandSideComputation const compute_velocity,
      Parameters const& parameters,
      Solution* solution) = 0;

};

}  // namespace integrators
}  // namespace principia

#include "integrators/integrator_body.hpp"