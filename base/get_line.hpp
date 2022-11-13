#pragma once

#include <fstream>
#include <string>

#include "base/not_null.hpp"

namespace principia {
namespace base {
namespace internal_get_line {

// Recursively reads a line of arbitrary length.
std::string GetLine(std::ifstream& stream);

}  // namespace internal_get_line

using internal_get_line::GetLine;

}  // namespace base
}  // namespace principia

#include "base/get_line_body.hpp"
