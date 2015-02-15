#pragma once

#include <memory>

#include "base/not_null.hpp"
#include "physics/n_body_system.hpp"
#include "testing_utilities/solar_system.hpp"

namespace principia {

using base::not_null;
using physics::NBodySystem;
using testing_utilities::SolarSystem;

namespace benchmarks {

// Simulates the given |solar_system| for 100 years with a 45 min time step.
void SimulateSolarSystem(not_null<SolarSystem*> const solar_system);

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/n_body_system_body.hpp"
