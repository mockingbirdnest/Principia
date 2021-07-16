
#pragma once

#include "base/macros.hpp"

namespace principia {
namespace base {

extern char const BuildDate[];
extern char const Version[];

#pragma message("Compiler version: " STRINGIFY_EXPANSION(_MSC_FULL_VER))

}  // namespace base
}  // namespace principia
