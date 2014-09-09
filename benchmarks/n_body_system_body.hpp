#pragma once

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "physics/n_body_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system.hpp"


using principia::astronomy::JulianYear;
using principia::integrators::SPRKIntegrator;
using principia::physics::NBodySystem;
using principia::quantities::Length;
using principia::quantities::Speed;
using principia::si::Minute;
using principia::testing_utilities::ICRFJ2000Ecliptic;
using principia::testing_utilities::SolarSystem;

namespace principia {
namespace benchmarks {

void SimulateSolarSystem(SolarSystem* solar_system) {
  std::unique_ptr<NBodySystem<ICRFJ2000Ecliptic>> n_body_system(
      new NBodySystem<ICRFJ2000Ecliptic>());
  auto const trajectories = solar_system->trajectories();
  SPRKIntegrator<Length, Speed> integrator;
  integrator.Initialize(integrator.Order5Optimal());
  n_body_system->Integrate(integrator,
                           trajectories.front()->last_time() +
                               100 * JulianYear,              // t_max
                           45 * Minute,                       // Δt
                           0,                                 // sampling_period
                           trajectories);
}


}  // namespace benchmarks
}  // namespace principia
