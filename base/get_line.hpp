#pragma once

#include <fstream>
#include <string>

namespace principia {
namespace base {
namespace _get_line {
namespace internal {

// Recursively reads a line of arbitrary length.
std::string GetLine(std::ifstream& stream);

}  // namespace internal

using internal::GetLine;

}  // namespace _get_line
}  // namespace base
}  // namespace principia

#include "base/get_line_body.hpp"
