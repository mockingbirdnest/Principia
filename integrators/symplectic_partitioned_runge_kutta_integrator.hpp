#pragma once

#include <vector>

#include "integrators/symplectic_integrator.hpp"

namespace principia {
namespace integrators {

class SPRKIntegrator : public SymplecticIntegrator {
 public:
  SPRKIntegrator();
  ~SPRKIntegrator() override = default;

  std::vector<std::vector<double>> const& Order5Optimal() const;

  void Initialize(Coefficients const& coefficients) override;

  void Solve(RightHandSideComputation const compute_force,
             AutonomousRightHandSideComputation const compute_velocity,
             Parameters const& parameters,
             Solution* solution) override;

 private:
  // The tableau.
  int stages_;
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> c_;
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
