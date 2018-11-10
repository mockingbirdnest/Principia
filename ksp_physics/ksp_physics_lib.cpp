
#include "base/macros.hpp"

#if OS_WIN

#include "ksp_physics/ksp_physics_lib.hpp"

#include <windows.h>
#include <psapi.h>

#include "base/version.hpp"
#include "glog/logging.h"

namespace principia {
namespace physics {

void LogPhysicsDLLBaseAddress() {
  LOG(INFO) << "Principia KSP physics DLL version " << principia::base::Version
            << " built on " << principia::base::BuildDate
            << " by " << principia::base::CompilerName
            << " version " << principia::base::CompilerVersion
            << " for " << principia::base::OperatingSystem
            << " " << principia::base::Architecture;
  MODULEINFO module_info;
  memset(&module_info, 0, sizeof(module_info));
  CHECK(GetModuleInformation(GetCurrentProcess(),
                             GetModuleHandle(TEXT("physics")),
                             &module_info,
                             sizeof(module_info)));
  LOG(INFO) << "Base address is " << module_info.lpBaseOfDll;
}

}  // namespace physics
}  // namespace principia

#endif  // OS_WIN
