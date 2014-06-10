#pragma once

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "physics/n_body_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system.hpp"

namespace principia {
namespace benchmarks {

physics::NBodySystem<testing_utilities::ICRFJ2000EclipticFrame> *
SimulateSolarSystem() {
  physics::NBodySystem<testing_utilities::ICRFJ2000EclipticFrame>* system =
    testing_utilities::SolarSystemAtSputnikLaunch();
  integrators::SPRKIntegrator integrator;
  integrator.Initialize(integrator.Order5Optimal());
  system->Integrate(
      integrator,
      testing_utilities::SputnikLaunchDate +
          100 * astronomy::JulianYear,       // t_max
      45 * si::Minute,                       // Δt
      0);                                    // sampling_period
  return system;
}


}  // namespace benchmarks
}  // namespace principia
