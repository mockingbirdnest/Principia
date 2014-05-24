#pragma once

#include "integrators/integrator.hpp"

namespace principia {
namespace integrators {

class SPRKIntegrator : public Integrator {
 public:
  SPRKIntegrator();
  virtual ~SPRKIntegrator();

  virtual void Increment(
      RightHandSideComputation const compute_force,
      AutonomousRightHandSideComputation const compute_velocity,
      Parameters const& parameters,
      Solution* solution);

};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"