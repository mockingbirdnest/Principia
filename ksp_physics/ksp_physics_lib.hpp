#pragma once

#include "base/macros.hpp"
#include "ksp_plugin/frames.hpp"

// If we are exporting, load the definitions; otherwise, just declare an extern
// instantiation.
#if !PHYSICS_DLL_IMPORT
#include "physics/ephemeris.hpp"
#include "physics/continuous_trajectory.hpp"
#endif

namespace principia {
namespace physics {

#if OS_WIN
PHYSICS_DLL void LogPhysicsDLLBaseAddress();
#endif

namespace internal_ephemeris {

PHYSICS_DLL_TEMPLATE_CLASS ContinuousTrajectory<ksp_plugin::Barycentric>;
PHYSICS_DLL_TEMPLATE_CLASS Ephemeris<ksp_plugin::Barycentric>;

}  // namespace internal_ephemeris
}  // namespace physics
}  // namespace principia
