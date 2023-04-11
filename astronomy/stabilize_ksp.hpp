#pragma once

#include "physics/solar_system.hpp"

namespace principia {
namespace astronomy {
namespace _stabilize_ksp {
namespace internal {

using namespace principia::physics::_solar_system;

// Patches the given |solar_system|, which is expected to be the stock KSP, to
// make the Jool system stable.
template<typename Frame>
void StabilizeKSP(SolarSystem<Frame>& solar_system);

}  // namespace internal

using internal::StabilizeKSP;

}  // namespace _stabilize_ksp
}  // namespace astronomy
}  // namespace principia

#include "astronomy/stabilize_ksp_body.hpp"
