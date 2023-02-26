#pragma once

#include <fstream>
#include <string>

#include "base/not_null.hpp"

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

namespace principia::base {
using namespace principia::base::_get_line;
}  // namespace principia::base

#include "base/get_line_body.hpp"
