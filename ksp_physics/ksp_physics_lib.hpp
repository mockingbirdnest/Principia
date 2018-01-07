#pragma once

#include "base/macros.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/ephemeris.hpp"

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
