
#pragma once

#include <fstream>
#include <string>

#include "base/not_null.hpp"

namespace principia {
namespace base {

// Recursively reads a line of arbitrary length.
std::string GetLine(not_null<std::ifstream*> stream);

}  // namespace base
}  // namespace principia

#include "base/get_line_body.hpp"
