
#pragma once

#include <cstdlib>
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
char const* const CompilerName = "Clang-cl";
char const* const CompilerVersion = __VERSION__;
#elif defined(__clang__)
#define PRINCIPIA_COMPILER_CLANG 1
char const* const CompilerName = "Clang";
char const* const CompilerVersion = __VERSION__;
#elif defined(_MSC_VER)
#define PRINCIPIA_COMPILER_MSVC 1
char const* const CompilerName = "Microsoft Visual C++";
char const* const CompilerVersion = STRINGIFY_EXPANSION(_MSC_FULL_VER);
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#define PRINCIPIA_COMPILER_ICC 1
char const* const CompilerName = "Intel C++ Compiler";
char const* const CompilerVersion = __VERSION__;
#elif defined(__GNUC__)
#define PRINCIPIA_COMPILER_GCC 1
char const* const CompilerName = "G++";
char const* const CompilerVersion = __VERSION__;
#else
#error "What is this, Borland C++?"
#endif

#if defined(__APPLE__)
#define OS_MACOSX 1
char const* const OperatingSystem = "OS X";
#elif defined(__linux__)
#define OS_LINUX 1
char const* const OperatingSystem = "Linux";
#elif defined(__FreeBSD__)
#define OS_FREEBSD 1
char const* const OperatingSystem = "FreeBSD";
#elif defined(_WIN32)
#define OS_WIN 1
char const* const OperatingSystem = "Windows";
#else
#error "Try OS/360."
#endif

#if defined(__i386) || defined(_M_IX86)
#define ARCH_CPU_X86_FAMILY 1
#define ARCH_CPU_X86 1
#define ARCH_CPU_32_BITS 1
#define ARCH_CPU_LITTLE_ENDIAN 1
char const* const Architecture = "x86";
#elif defined(_M_X64) || defined(__x86_64__)
#define ARCH_CPU_X86_FAMILY 1
#define ARCH_CPU_X86_64 1
#define ARCH_CPU_64_BITS 1
#define ARCH_CPU_LITTLE_ENDIAN 1
char const* const Architecture = "x86-64";
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
#if defined(PRINCIPIA_DLL)
#  error "PRINCIPIA_DLL already defined"
#else
#  if OS_WIN
#    if PRINCIPIA_DLL_IMPORT
#      define PRINCIPIA_DLL __declspec(dllimport)
#    else
#      define PRINCIPIA_DLL __declspec(dllexport)
#    endif
#  else
#    define PRINCIPIA_DLL __attribute__((visibility("default")))
#  endif
#endif

// DLL-exported functions for isolated physics optimization.
#if defined(PHYSICS_DLL)
#  error "PHYSICS_DLL already defined"
#else
#  if OS_WIN
#    if PHYSICS_DLL_IMPORT
#      define PHYSICS_DLL __declspec(dllimport)
#      define PHYSICS_DLL_TEMPLATE_CLASS extern template class PHYSICS_DLL
#    else
#      define PHYSICS_DLL __declspec(dllexport)
#      define PHYSICS_DLL_TEMPLATE_CLASS template class PHYSICS_DLL
#    endif
#  else
#    define PHYSICS_DLL __attribute__((visibility("default")))
     // No isolated physics library on Linux or Macintosh at the moment.
#    define PHYSICS_DLL_TEMPLATE_CLASS template class PHYSICS_DLL
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
inline void noreturn() { std::exit(0); }

// Used to force inlining.
#if PRINCIPIA_COMPILER_CLANG    ||  \
    PRINCIPIA_COMPILER_CLANG_CL ||  \
    PRINCIPIA_COMPILER_GCC
#  define FORCE_INLINE(specifiers) [[gnu::always_inline]] specifiers // NOLINT
#elif PRINCIPIA_COMPILER_MSVC
#  define FORCE_INLINE(specifiers) specifiers __forceinline
#elif PRINCIPIA_COMPILER_ICC
#  define FORCE_INLINE(specifiers) __attribute__((always_inline))
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

// We assume that the processor is at least a Prescott since we only support
// 64-bit architectures.
#define PRINCIPIA_USE_SSE3_INTRINSICS !_DEBUG

// Thread-safety analysis.
#if PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL
#  define THREAD_ANNOTATION_ATTRIBUTE__(x) __attribute__((x))
#  define EXCLUDES(...) \
       THREAD_ANNOTATION_ATTRIBUTE__(locks_excluded(__VA_ARGS__))
#  define GUARDED_BY(...) \
       THREAD_ANNOTATION_ATTRIBUTE__(guarded_by(__VA_ARGS__))
#  define REQUIRES(...) \
       THREAD_ANNOTATION_ATTRIBUTE__(requires_capability(__VA_ARGS__))
#  define REQUIRES_SHARED(...) \
       THREAD_ANNOTATION_ATTRIBUTE__(requires_shared_capability(__VA_ARGS__))
#else
#  define EXCLUDES(x)
#  define GUARDED_BY(x)
#  define REQUIRES(x)
#  define REQUIRES_SHARED(x)
#endif

// Unicode.
#if OS_WIN
#  define UNICODE_PATH(x) u ## x
#else
#  define UNICODE_PATH(x) u8 ## x
#endif

// Mutex.
#if !OS_MACOSX
#  define HAS_SHARED_MUTEX 1
#endif

#define NAMED(expression) #expression << ": " << (expression)

// A macro to allow glog checking within C++11 constexpr code.  If |condition|
// is true, evaluates to |expression|.  Otherwise, results in a CHECK failure at
// runtime and a compilation error due to a call to non-constexpr code at
// compile time.
// NOTE(egg): in the failure case, the |LOG(FATAL)| is wrapped in a lambda.  The
// reason is that |LOG(FATAL)| constructs a |google::LogMessageFatal|, and
// |google::LogMessage::Fail()| is called in its destructor.  As a temporary,
// the |google::LogMessageFatal| is destroyed as the last step in evaluating the
// enclosing full-expression.  If we simply used the comma operator, the entire
// ternary |((condition) ? (expression) : (CHECK(condition), (expression)))|
// would be part of the enclosing full-expression, so that |expression| would
// get evaluated before the |CHECK| failure, possibly triggering all sorts of
// terrible UB or other checks (|DateDeathTest| provides a couple of examples).
// With the lambda, the full-expression forms the expression statement
// |LOG(FATAL) << "Check failed: " #condition " ";|, so that failure occurs
// before we return from the lambda's function call operator, and |expression|
// is never evaluated.  We do not use |CHECK| because that would require
// capture, but this should produce the same output.
#define CHECKING(condition, expression)                                      \
  ((condition) ? (expression)                                                \
               : (([] { LOG(FATAL) << "Check failed: " #condition " "; })(), \
                  (expression)))

// Clang for some reason doesn't like FP arithmetic that yields infinities in
// constexpr code (MSVC and GCC are fine with that).
#if PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL
#  define CONSTEXPR_INFINITY const
#else
#  define CONSTEXPR_INFINITY constexpr
#endif

#if PRINCIPIA_COMPILER_MSVC
#define MSVC_ONLY_TEST(test_name) test_name
#else
#define MSVC_ONLY_TEST(test_name) DISABLED_##test_name
#endif

// For templates in macro parameters.
#define TEMPLATE(...) template<__VA_ARGS__>

// Forward declaration of a class or struct declared in an internal namespace
// according to #602.
// Usage:
// FORWARD_DECLARE_FROM(p1, struct, T);
// FORWARD_DECLARE_FROM(p2, TEMPLATE(int i) class, U);
#define FORWARD_DECLARE_FROM(package_name,           \
                             template_and_class_key, \
                             declared_name)          \
namespace internal_##package_name {                  \
template_and_class_key declared_name;                \
}                                                    \
using internal_##package_name::declared_name

}  // namespace base
}  // namespace principia
