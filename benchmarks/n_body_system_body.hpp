#pragma once

#include "base/not_null.hpp"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "physics/n_body_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system.hpp"

namespace principia {

using astronomy::JulianYear;
using base::not_null;
using integrators::SPRKIntegrator;
using physics::NBodySystem;
using quantities::Length;
using quantities::Speed;
using si::Minute;
using testing_utilities::ICRFJ2000Ecliptic;
using testing_utilities::SolarSystem;

namespace benchmarks {

void SimulateSolarSystem(not_null<SolarSystem*> const solar_system) {
  auto const n_body_system = std::make_unique<NBodySystem<ICRFJ2000Ecliptic>>();
  auto const trajectories = solar_system->trajectories();
  SPRKIntegrator<Length, Speed> integrator;
  integrator.Initialize(integrator.McLachlanAtela1992Order5Optimal());
  n_body_system->Integrate(integrator,
                           trajectories.front()->last().time() +
                               100 * JulianYear,              // t_max
                           45 * Minute,                       // Δt
                           0,                                 // sampling_period
                           false,                             // tmax_is_exact
                           trajectories);
}

}  // namespace benchmarks
}  // namespace principia
