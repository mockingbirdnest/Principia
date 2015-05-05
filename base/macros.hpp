#pragma once

#include <string>

namespace principia {
namespace base {

#define STRINGIFY(X) #X
#define STRINGIFY_EXPANSION(X) STRINGIFY(X)

// See http://goo.gl/2EVxN4 for a partial overview of compiler detection and
// version macros.  We cannot use |COMPILER_MSVC| because it conflicts with
// a macro in the benchmark library, so the macros have obnoxiously long names.
// TODO(phl): See whether that |COMPILER_MSVC| macro can be removed from port.h.
#if defined(_MSC_VER) && defined(__clang__)
#define PRINCIPIA_COMPILER_CLANG_CL 1
char const* const kCompilerName = "Clang-cl";
char const* const kCompilerVersion = __VERSION__;
#elif defined(__clang__)
#define PRINCIPIA_COMPILER_CLANG 1
char const* const kCompilerName = "Clang";
char const* const kCompilerVersion = __VERSION__;
#elif defined(_MSC_VER)
#define PRINCIPIA_COMPILER_MSVC 1
char const* const kCompilerName = "Microsoft Visual C++";
char const* const kCompilerVersion = STRINGIFY_EXPANSION(_MSC_FULL_VER);
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#define PRINCIPIA_COMPILER_ICC 1
char const* const kCompilerName = "Intel C++ Compiler";
char const* const kCompilerVersion = __VERSION__;
#elif defined(__GNUC__)
#define PRINCIPIA_COMPILER_GCC 1
char const* const kCompilerName = "G++";
char const* const kCompilerVersion = __VERSION__;
#else
#error "What is this, Borland C++?"
#endif

#if defined(__APPLE__)
#define OS_MACOSX 1
char const* const kOperatingSystem = "OS X";
#elif defined(__linux__)
#define OS_LINUX 1
char const* const kOperatingSystem = "Linux";
#elif defined(__FreeBSD__)
#define OS_FREEBSD 1
char const* const kOperatingSystem = "FreeBSD";
#elif defined(_WIN32)
#define OS_WIN 1
char const* const kOperatingSystem = "Windows";
#else
#error "Try OS/360."
#endif

#if defined(__i386) || defined(_M_IX86)
#define ARCH_CPU_X86_FAMILY 1
#define ARCH_CPU_X86 1
#define ARCH_CPU_32_BITS 1
#define ARCH_CPU_LITTLE_ENDIAN 1
char const* const kArchitecture = "x86";
#elif defined(_M_X64) || defined(__x86_64__)
#define ARCH_CPU_X86_FAMILY 1
#define ARCH_CPU_X86_64 1
#define ARCH_CPU_64_BITS 1
#define ARCH_CPU_LITTLE_ENDIAN 1
char const* const kArchitecture = "x86-64";
#else
#error "Have you tried a Cray-1?"
#endif

#if defined(CDECL)
#  error "CDECL already defined"
#else
// Architecture macros from http://goo.gl/ZypnO8.
// We use cdecl on x86, the calling convention is unambiguous on x86-64.
#  if ARCH_CPU_X86
#    if PRINCIPIA_COMPILER_CLANG ||  \
        PRINCIPIA_COMPILER_MSVC ||   \
        PRINCIPIA_COMPILER_CLANG_CL
#      define CDECL __cdecl
#    elif PRINCIPIA_COMPILER_ICC || PRINCIPIA_COMPILER_GCC
#      define CDECL __attribute__((cdecl))
#    else
#      error "Get a real compiler!"
#    endif
#  elif ARCH_CPU_X86_64
#    define CDECL
#  else
#    error "Have you tried a Cray-1?"
#  endif
#endif

// DLL-exported functions for interfacing with Platform Invocation Services.
#if defined(DLLEXPORT)
#  error "DLLEXPORT already defined"
#else
#  if OS_WIN
#    define DLLEXPORT __declspec(dllexport)
#  else
#    define DLLEXPORT __attribute__((visibility("default")))
#  endif
#endif

// A function for use on control paths that don't return a value, typically
// because they end with a |LOG(FATAL)|.
#if PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL
[[noreturn]]
#elif PRINCIPIA_COMPILER_MSVC
__declspec(noreturn)
#elif PRINCIPIA_COMPILER_ICC
__attribute__((noreturn))
#else
#error "What compiler is this?"
#endif
inline void noreturn() { exit(0); }

// Used to force inlining.
#if PRINCIPIA_COMPILER_CLANG    ||  \
    PRINCIPIA_COMPILER_CLANG_CL ||  \
    PRINCIPIA_COMPILER_GCC
#  define FORCE_INLINE [[gnu::always_inline]]  // NOLINT(whitespace/braces)
#elif PRINCIPIA_COMPILER_MSVC
#  define FORCE_INLINE __forceinline
#elif PRINCIPIA_COMPILER_ICC
#  define FORCE_INLINE __attribute__((always_inline))
#else
#  error "What compiler is this?"
#endif

// Used to emit the function signature.
#if PRINCIPIA_COMPILER_CLANG    ||  \
    PRINCIPIA_COMPILER_CLANG_CL ||  \
    PRINCIPIA_COMPILER_GCC
#  define FUNCTION_SIGNATURE __PRETTY_FUNCTION__
#elif PRINCIPIA_COMPILER_MSVC
#  define FUNCTION_SIGNATURE __FUNCSIG__
#else
#  error "What compiler is this?"
#endif

// Thread-safety analysis.
#if PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL
#  define THREAD_ANNOTATION_ATTRIBUTE__(x) __attribute__((x))
#  define GUARDED_BY(x) THREAD_ANNOTATION_ATTRIBUTE__(guarded_by(x))
#else
#  define GUARDED_BY(x)
#endif

// A workaround for a MSVC bug wherein a |typename| is required by the standard
// and by clang but forbidden by MSVC.
#if PRINCIPIA_COMPILER_MSVC
#  define TYPENAME
#else
#  define TYPENAME typename
#endif

// Same as above, with |template|.
#if PRINCIPIA_COMPILER_MSVC
#  define TEMPLATE
#else
#  define TEMPLATE template
#endif

#define VLOG_AND_RETURN(verboselevel, expression)                  \
  do {                                                             \
    auto const& value__ = (expression);                            \
    VLOG(verboselevel) << __FUNCTION__ << " returns " << value__;  \
    return value__;                                                \
  } while (false)

#define NAMED(expression) #expression << ": " << (expression)

}  // namespace base
}  // namespace principia
