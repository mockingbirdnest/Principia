#pragma once

#include "physics/solar_system.hpp"

namespace principia {
namespace astronomy {
namespace ksp_stabilizer_internal {

using physics::SolarSystem;

// Patches the given |solar_system|, which is expected to be the stock KSP, to
// make the Jool system stable.
template<typename Frame>
void KSPStabilizer(SolarSystem<Frame>& solar_system);

}  // namespace ksp_stabilizer_internal

using ksp_stabilizer_internal::KSPStabilizer;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/ksp_stabilizer_body.hpp"
