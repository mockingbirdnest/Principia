#pragma once

// DLL-exported functions for interfacing with Platform Invocation Services.

#if defined(_WIN32)
#define DLLEXPORT __declspec(dllexport)
#elif defined(__INTEL_COMPILER) || defined(__GNUC__) || defined(__clang__)
#define DLLEXPORT __attribute__((visibility("default")))
#else
#error "What compiler do you think you're using?"
#endif

namespace principia {
namespace ksp_plugin {

extern "C" DLLEXPORT
int Say33();

extern "C"
__declspec(dllexport)
char const* SayHello();

}  // namespace ksp_plugin
}  // namespace principia
