#pragma once

#include "integrators/integrator.hpp"

namespace principia {
namespace integrators {

class SPRKIntegrator : public Integrator {
 public:
  static std::vector<std::vector<double>> const order_5_optimal = {
      { 0.339839625839110000,
       -0.088601336903027329,
        0.5858564768259621188,
       -0.603039356536491888,
        0.3235807965546976394,
        0.4423637942197494587},
      { 0.1193900292875672758,
        0.6989273703824752308,
       -0.1713123582716007754,
        0.4012695022513534480,
        0.0107050818482359840,
       -0.0589796254980311632}};

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