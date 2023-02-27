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
# if _HAS_CXX20
# define PRINCIPIA_COMPILER_MSVC_HAS_CXX20 1
# endif
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
#define PRINCIPIA_USE_FMA_IF_AVAILABLE !_DEBUG

// Set this to 1 to test analytical series based on piecewise Poisson series.
#define PRINCIPIA_CONTINUOUS_TRAJECTORY_SUPPORTS_PIECEWISE_POISSON_SERIES 0

// Thread-safety analysis.
#if PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL
#  define THREAD_ANNOTATION_ATTRIBUTE__(x) __attribute__((x))
#  define EXCLUDES(...) \
       THREAD_ANNOTATION_ATTRIBUTE__(locks_excluded(__VA_ARGS__))
#  define REQUIRES(...) \
       THREAD_ANNOTATION_ATTRIBUTE__(requires_capability(__VA_ARGS__))
#  define REQUIRES_SHARED(...) \
       THREAD_ANNOTATION_ATTRIBUTE__(requires_shared_capability(__VA_ARGS__))
#else
#  define EXCLUDES(x)
#  define REQUIRES(x)
#  define REQUIRES_SHARED(x)
#endif

// Unicode.
#if OS_WIN
#  define PRINCIPIA_UNICODE_PATH(x) u ## x
#else
#  define PRINCIPIA_UNICODE_PATH(x) u8 ## x
#endif

#define NAMED(expression) u8 ## #expression << ": " << (expression)

// Needed to circumvent lint warnings in constexpr functions where CHECK_LT and
// friends cannot be used.
#define CONSTEXPR_CHECK(condition) CHECK(condition)
#define CONSTEXPR_DCHECK(condition) DCHECK(condition)

// Lexicographic comparison (v1, v2, v3) ≥ (w1, w2, w3).
#define VERSION_GE(v1, v2, v3, w1, w2, w3)         \
  ((v1) > (w1) || ((v1) == (w1) && (v2) > (w2)) || \
    ((v1) == (w1) && (v2) == (w2) && (v3) >= (w3)))

#define CLANG_VERSION_GE(major, minor, patchlevel) \
  VERSION_GE(__clang_major__,                      \
             __clang_minor__,                      \
             __clang_patchlevel__,                 \
             major,                                \
             minor,                                \
             patchlevel)

// Clang does not like FP arithmetic that yields a NaN in constexpr code.
// https://github.com/llvm/llvm-project/blob/llvmorg-13.0.0/clang/lib/AST/ExprConstant.cpp#L2860-L2867
#if PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL
#  define CONSTEXPR_NAN const
#else
#  define CONSTEXPR_NAN constexpr
#endif

#if PRINCIPIA_COMPILER_MSVC
#define MSVC_ONLY_TEST(test_name) test_name
#else
#define MSVC_ONLY_TEST(test_name) DISABLED_##test_name
#endif

// For templates in macro parameters.
#define TEMPLATE(...) template<__VA_ARGS__>

// For circumventing
// https://developercommunity.visualstudio.com/content/problem/1256363/operator-call-incorrectly-marked-as-ambiguous-with.html.
#if PRINCIPIA_COMPILER_MSVC_HAS_CXX20
#define PRINCIPIA_MAX(l, r) ((l) > (r) ? (l) : (r))
#define PRINCIPIA_MAX3(x1, x2, x3) \
  PRINCIPIA_MAX((x1), PRINCIPIA_MAX((x2), (x3)))
#define PRINCIPIA_MAX4(x1, x2, x3, x4) \
  PRINCIPIA_MAX((x1), PRINCIPIA_MAX((x2), PRINCIPIA_MAX((x3), (x4))))
#else
#define PRINCIPIA_MAX(l, r) std::max((l), (r))
#define PRINCIPIA_MAX3(x1, x2, x3) std::max({(x1), (x2), (x3)})
#define PRINCIPIA_MAX4(x1, x2, x3, x4) std::max({(x1), (x2), (x3), (x4)})
#endif

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

#define FORWARD_DECLARE_FROM_NEW(                        \
    package_name, template_and_class_key, declared_name) \
  namespace _##package_name {                            \
    namespace internal {                                 \
    template_and_class_key declared_name;                \
    }                                                    \
    using internal::declared_name;                       \
  }

#define FORWARD_DECLARE_FUNCTION_FROM(package_name,        \
                                      template_and_result, \
                                      declared_name,       \
                                      parameters)          \
namespace internal_##package_name {                        \
template_and_result declared_name parameters;              \
}                                                          \
using internal_##package_name::declared_name

}  // namespace base
}  // namespace principia
