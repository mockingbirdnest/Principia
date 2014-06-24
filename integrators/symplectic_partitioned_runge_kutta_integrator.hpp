#pragma once

#include <vector>

#include "integrators/symplectic_integrator.hpp"

namespace principia {
namespace integrators {

template<typename Position, typename Momentum>
class SPRKIntegrator : public SymplecticIntegrator<Position, Momentum> {
 public:
  typedef typename SymplecticIntegrator<Position, Momentum>::Coefficients
      Coefficients;
  typedef typename SymplecticIntegrator<Position, Momentum>::Parameters
      Parameters;
  typedef typename SymplecticIntegrator<Position, Momentum>::SystemState
      SystemState;

  SPRKIntegrator();
  ~SPRKIntegrator() override = default;

  std::vector<std::vector<double>> const& Order5Optimal() const;

  void Initialize(Coefficients const& coefficients) override;

  template<typename AutonomousRightHandSideComputation,
           typename RightHandSideComputation>
  void Solve(RightHandSideComputation compute_force,
             AutonomousRightHandSideComputation compute_velocity,
             Parameters const& parameters,
             std::vector<SystemState>* solution) const;

 private:
  int stages_;

  // The position and momentum nodes.
  std::vector<double> a_;
  std::vector<double> b_;

  // The weights.
  std::vector<double> c_;
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
