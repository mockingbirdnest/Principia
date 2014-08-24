#pragma once

#include <memory>

#include "physics/n_body_system.hpp"
#include "testing_utilities/solar_system.hpp"

using principia::physics::NBodySystem;
using principia::testing_utilities::SolarSystem;

namespace principia {
namespace benchmarks {

// Simulates the given |solar_system| for 100 years with a 45 min time step.
void SimulateSolarSystem(SolarSystem* solar_system);

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/n_body_system_body.hpp"
