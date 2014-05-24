#pragma once

#include <vector>

#include "integrators/symplectic_integrator.hpp"

namespace principia {
namespace integrators {

class SPRKIntegrator : public SymplecticIntegrator {
 public:
  SPRKIntegrator() = default;
  ~SPRKIntegrator() override = default;

  std::vector<std::vector<double>> const& Order5Optimal() const;

  void Solve(RightHandSideComputation const compute_force,
             AutonomousRightHandSideComputation const compute_velocity,
             Parameters const& parameters,
             Solution* solution) override;

};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
