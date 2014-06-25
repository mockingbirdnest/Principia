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
using principia::testing_utilities::kСпутникLaunchDate;
using principia::testing_utilities::SolarSystemAtСпутникLaunch;

namespace principia {
namespace benchmarks {

std::unique_ptr<NBodySystem<ICRFJ2000EclipticFrame>> SimulateSolarSystem() {
  std::unique_ptr<NBodySystem<ICRFJ2000EclipticFrame>> system =
      SolarSystemAtСпутникLaunch();
  SPRKIntegrator<Length, Speed> integrator;
  integrator.Initialize(integrator.Order5Optimal());
  system->Integrate(integrator,
                    kСпутникLaunchDate + 100 * JulianYear,  // t_max
                    45 * Minute,                           // Δt
                    0);                                    // sampling_period
  return system;
}


}  // namespace benchmarks
}  // namespace principia
