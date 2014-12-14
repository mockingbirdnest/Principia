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

// A function for use on control paths that don't return a value, typically
// because they end with a |LOG(FATAL)|.
#ifdef __clang__
[[noreturn]]
#elif _MSC_VER
__declspec(noreturn)
#endif
inline void noreturn() { exit(0); }

// Used to force inlining.
#ifdef __clang__
#define FORCE_INLINE [[gnu::always_inline]]  // NOLINT(whitespace/braces)
#elif _MSC_VER
#define FORCE_INLINE __forceinline
#endif

#define VLOG_AND_RETURN(verboselevel, expression)                  \
  do {                                                             \
    auto const& value__ = (expression);                            \
    VLOG(verboselevel) << __FUNCTION__ << " returns " << value__;  \
    return value__;                                                \
  } while (false)

}  // namespace base
}  // namespace principia
