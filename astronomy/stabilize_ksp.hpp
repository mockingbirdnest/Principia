#pragma once

#include "physics/solar_system.hpp"

namespace principia {
namespace astronomy {
namespace stabilize_ksp_internal {
namespace _stabilize_ksp {
namespace internal {

using physics::SolarSystem;

// Patches the given |solar_system|, which is expected to be the stock KSP, to
// make the Jool system stable.
template<typename Frame>
void StabilizeKSP(SolarSystem<Frame>& solar_system);

}  // namespace internal

using internal::SolarSystem;
using internal::StabilizeKSP;

}  // namespace _stabilize_ksp
}  // namespace stabilize_ksp_internal

using stabilize_ksp_internal::StabilizeKSP;

}  // namespace astronomy
}  // namespace principia

namespace principia::astronomy {
using namespace principia::astronomy::_stabilize_ksp;
}  // namespace principia::astronomy

#include "astronomy/stabilize_ksp_body.hpp"
