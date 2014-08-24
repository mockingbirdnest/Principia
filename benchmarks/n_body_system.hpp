#pragma once

#include <memory>

#include "physics/n_body_system.hpp"
#include "testing_utilities/solar_system.hpp"

using principia::physics::NBodySystem;
using principia::testing_utilities::SolarSystem;

namespace principia {
namespace benchmarks {

// Simulates the solar system as at Спутник launch for 100 years with a 45 min
// time step.  The caller gets ownership of the returned object.  The |state|
// object may be used to tune benchmarking.
std::unique_ptr<SolarSystem> SimulateSolarSystem(benchmark::State* state);

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/n_body_system_body.hpp"
