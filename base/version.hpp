#pragma once

#include "base/macros.hpp"

namespace principia {
namespace base {

extern char const BuildDate[];
extern char const Version[];

#if OS_WIN
#pragma message("Compiler version: " STRINGIFY_EXPANSION(_MSC_FULL_VER))
#endif

}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_version;
}  // namespace principia::base
