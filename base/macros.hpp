#pragma once

#include <cstdlib>
#include <string>

namespace principia {
namespace base {

#define STRINGIFY(X) #X
#define STRINGIFY_EXPANSION(X) STRINGIFY(X)

#define PRINCIPIA_CONCATENATE_SENTINEL(X) X##0x4C##45##4E##49##54##4E##45##53
// True if X is #defined to nothing, false if X is #defined to an identifier or
// is not defined.  This macro should not be used with macros that expand to
// something other than an identifier.
#define PRINCIPIA_MACRO_IS_EMPTY(X)                                     \
    (PRINCIPIA_CONCATENATE_SENTINEL(X) ==                               \
     ('S' << 000 | 'E' << 010 | 'N' << 020 | 'T' << 030 | 'I' << 040 |  \
      'N' << 050 | 'E' << 060 | 'L' << 070))

// See http://goo.gl/2EVxN4 for a partial overview of compiler detection and
// version macros.
#if defined(_MSC_VER) && defined(__clang__)
#define PRINCIPIA_COMPILER_CLANG_CL 1
#define PRINCIPIA_VECTORCALL
char const* const CompilerName = "Clang-cl";
char const* const CompilerVersion = __VERSION__;
#elif defined(__clang__)
#define PRINCIPIA_COMPILER_CLANG 1
#define PRINCIPIA_VECTORCALL
char const* const CompilerName = "Clang";
char const* const CompilerVersion = __VERSION__;
#elif defined(_MSC_VER)
#define PRINCIPIA_VECTORCALL __vectorcall
#define PRINCIPIA_COMPILER_MSVC 1
char const* const CompilerName = "Microsoft Visual C++";
char const* const CompilerVersion = STRINGIFY_EXPANSION(_MSC_FULL_VER);
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#define PRINCIPIA_COMPILER_ICC 1
#define PRINCIPIA_VECTORCALL
char const* const CompilerName = "Intel C++ Compiler";
char const* const CompilerVersion = __VERSION__;
#elif defined(__GNUC__)
#define PRINCIPIA_COMPILER_GCC 1
#define PRINCIPIA_VECTORCALL
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
#define PRINCIPIA_USE_SSE3_INTRINSICS() !_DEBUG
#define PRINCIPIA_USE_FMA_IF_AVAILABLE() !_DEBUG

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

#define PRINCIPIA_REPEAT8(...) \
  PRINCIPIA_REPEAT4(__VA_ARGS__) PRINCIPIA_REPEAT4(__VA_ARGS__)
#define PRINCIPIA_REPEAT4(...) \
  PRINCIPIA_REPEAT2(__VA_ARGS__) PRINCIPIA_REPEAT2(__VA_ARGS__)
#define PRINCIPIA_REPEAT2(...) \
  PRINCIPIA_REPEAT1(__VA_ARGS__) PRINCIPIA_REPEAT1(__VA_ARGS__)
#define PRINCIPIA_REPEAT1(...) __VA_ARGS__

// The macro magic is inspired from http://jhnet.co.uk/articles/cpp_magic.  Note
// that we are using __VA_OPT__ to stop the recursion and detect empty argument
// lists because we are modern.
#define PRINCIPIA_EMPTY()
#define PRINCIPIA_DEFER1(m) m PRINCIPIA_EMPTY()

#define PRINCIPIA_EVAL(...) PRINCIPIA_EVAL16(__VA_ARGS__)
#define PRINCIPIA_EVAL16(...) PRINCIPIA_EVAL8(PRINCIPIA_EVAL8(__VA_ARGS__))
#define PRINCIPIA_EVAL8(...) PRINCIPIA_EVAL4(PRINCIPIA_EVAL4(__VA_ARGS__))
#define PRINCIPIA_EVAL4(...) PRINCIPIA_EVAL2(PRINCIPIA_EVAL2(__VA_ARGS__))
#define PRINCIPIA_EVAL2(...) PRINCIPIA_EVAL1(PRINCIPIA_EVAL1(__VA_ARGS__))
#define PRINCIPIA_EVAL1(...) __VA_ARGS__

// Applies m_LAST to its last argument and m_NOT_LAST to the others.
#define PRINCIPIA_MAP1_LAST(m, a0, ...) \
  m##__VA_OPT__(_NOT)##_LAST(a0)       \
      __VA_OPT__(PRINCIPIA_DEFER1(_PRINCIPIA_MAP1_LAST)()(m, __VA_ARGS__))
#define _PRINCIPIA_MAP1_LAST() PRINCIPIA_MAP1_LAST

#define PRINCIPIA_MAP2(m, a0, a1, ...) \
  m(a0, a1) __VA_OPT__(PRINCIPIA_DEFER1(_PRINCIPIA_MAP2)()(m, a0, __VA_ARGS__))
#define _PRINCIPIA_MAP2() PRINCIPIA_MAP2

#define USING_DIRECTIVE_INTO(from_package_name, into_package_name) \
  namespace into_package_name {                                    \
    namespace internal {                                           \
      using namespace from_package_name;                           \
    }                                                              \
  }

#define USING_DIRECTIVES_INTO(from_package_name, ...) \
  __VA_OPT__(PRINCIPIA_EVAL16(                        \
      PRINCIPIA_MAP2(USING_DIRECTIVE_INTO, from_package_name, __VA_ARGS__)))

// Forward declaration of a class or struct declared in an internal namespace
// according to #602.
// Usage:
//   FORWARD_DECLARE(struct, T, FROM(pack));
//   FORWARD_DECLARE(TEMPLATE(int i) class, U, FROM(pack));
//   FORWARD_DECLARE(class, T,
//                   FROM(pack_a, pack_b),
//                   INTO(pack1_a, pack1_b),
//                   INTO(pack2));
// There must be a FROM argument.  Optionally, there can be multiple INTO
// arguments, in which case a using directive for pack is inserted into pack1,
// pack2, etc.

#define NAMESPACE_NOT_LAST(x) x ::
#define NAMESPACE_LAST(x) _ ## x

#define FROM(...) PRINCIPIA_EVAL16(PRINCIPIA_MAP1_LAST(NAMESPACE, __VA_ARGS__))
#define INTO(...) PRINCIPIA_EVAL16(PRINCIPIA_MAP1_LAST(NAMESPACE, __VA_ARGS__))

#define FORWARD_DECLARE(                                           \
    template_and_class_key, declared_name, from_package_name, ...) \
  namespace from_package_name {                                    \
    namespace internal {                                           \
      template_and_class_key declared_name;                        \
    }                                                              \
    using internal::declared_name;                                 \
  }                                                                \
  USING_DIRECTIVES_INTO(from_package_name, __VA_ARGS__)

#define FORWARD_DECLARE_FUNCTION(                                           \
    template_and_result, declared_name, parameters, from_package_name, ...) \
  namespace from_package_name {                                             \
    namespace internal {                                                    \
      template_and_result declared_name parameters;                         \
    }                                                                       \
    using internal::declared_name;                                          \
  }                                                                         \
  USING_DIRECTIVES_INTO(from_package_name, __VA_ARGS__)

}  // namespace base
}  // namespace principia
