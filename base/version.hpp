#pragma once

#include "base/macros.hpp"  // 🧙 For OS_WIN.

namespace principia {
namespace base {
namespace _version {
namespace internal {

extern char const BuildDate[];
extern char const Version[];

#if PRINCIPIA_COMPILER_MSVC
#pragma message("Compiler version: " STRINGIFY_EXPANSION(_MSC_FULL_VER))
#elif PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL
#pragma message("Compiler version: " STRINGIFY_EXPANSION(__clang_major__) "." STRINGIFY_EXPANSION(__clang_minor__) "." STRINGIFY_EXPANSION(__clang_patchlevel__))  // NOLINT
#endif

}  // namespace internal

using internal::BuildDate;
using internal::Version;

}  // namespace _version
}  // namespace base
}  // namespace principia
