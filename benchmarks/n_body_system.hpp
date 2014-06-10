#pragma once

#include "physics/n_body_system.hpp"
#include "testing_utilities/solar_system.hpp"

namespace principia {
namespace benchmarks {

physics::NBodySystem<testing_utilities::ICRFJ2000EclipticFrame> *
SimulateSolarSystem();

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/n_body_system_body.hpp"
