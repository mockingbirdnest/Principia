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
using principia::testing_utilities::ICRFJ2000EclipticFrame;
using principia::testing_utilities::SolarSystem;

namespace principia {
namespace benchmarks {

std::unique_ptr<SolarSystem> SimulateSolarSystem(benchmark::State* state) {
  state->PauseTiming();
  std::unique_ptr<SolarSystem> solar_system(new SolarSystem);
  std::unique_ptr<NBodySystem<ICRFJ2000EclipticFrame>> n_body_system(
      new NBodySystem<ICRFJ2000EclipticFrame>(
          solar_system->massive_bodies(),
          solar_system->massless_bodies()));
  state->ResumeTiming();
  SPRKIntegrator<Length, Speed> integrator;
  integrator.Initialize(integrator.Order5Optimal());
  n_body_system->Integrate(integrator,
                           solar_system->спутник_launch_time() +
                               100 * JulianYear,              // t_max
                           45 * Minute,                       // Δt
                           0,                                 // sampling_period
                           solar_system->trajectories_at_спутник_launch_time());
  return solar_system;
}


}  // namespace benchmarks
}  // namespace principia
