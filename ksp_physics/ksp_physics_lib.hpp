#pragma once

#include "base/macros.hpp"
#include "ksp_plugin/frames.hpp"

// If we are exporting, load the definitions; otherwise, just declare an extern
// specialization.
#if !PRINCIPIA_DLL_IMPORT
#include "physics/ephemeris.hpp"
#endif

namespace principia {
namespace physics {

#if OS_WIN
PHYSICS_DLL void LogPhysicsDLLBaseAddress();
#endif

namespace internal_ephemeris {

PHYSICS_DLL_TEMPLATE_CLASS Ephemeris<ksp_plugin::Barycentric>;

}  // namespace internal_ephemeris
}  // namespace physics
}  // namespace principia
