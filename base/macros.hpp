#pragma once

namespace principia {
namespace base {

#if defined(CDECL)
#  error "CDECL already defined"
#else
// Architecture macros from http://goo.gl/ZypnO8.
// We use cdecl on x86, the calling convention is unambiguous on x86-64.
#  if defined(__i386) || defined(_M_IX86)
#    if defined(_MSC_VER) || defined(__clang__)
#      define CDECL __cdecl
#    elif defined(__GNUC__) || defined(__INTEL_COMPILER)
#      define CDECL __attribute__((cdecl))
#    else
#      error "Get a real compiler"
#    endif
#  elif defined(_M_X64) || defined(__x86_64__)
#    define CDECL
#  else
#    error "Have you tried a Cray-1?"
#  endif
#endif

// DLL-exported functions for interfacing with Platform Invocation Services.
#if defined(DLLEXPORT)
#  error "DLLEXPORT already defined"
#else
#  if defined(_WIN32) || defined(_WIN64)
#    define DLLEXPORT __declspec(dllexport)
#  else
#    define DLLEXPORT __attribute__((visibility("default")))
#  endif
#endif

}  // namespace base
}  // namespace principia