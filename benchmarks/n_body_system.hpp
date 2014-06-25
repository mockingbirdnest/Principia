#pragma once

#include <memory>

#include "physics/n_body_system.hpp"
#include "testing_utilities/solar_system.hpp"

using principia::physics::NBodySystem;
using principia::testing_utilities::ICRFJ2000EclipticFrame;

namespace principia {
namespace benchmarks {

// Simulates the solar system as obtained by |SolarSystemAtSputnikLaunch| for
// 100 years with a 45 min time step.
// The caller gets ownership of the returned object.
std::unique_ptr<NBodySystem<ICRFJ2000EclipticFrame>> SimulateSolarSystem();

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/n_body_system_body.hpp"
