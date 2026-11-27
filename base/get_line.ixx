export module principia.base.get_line;

import std;

import principia.base.get_line;

namespace principia {
namespace base {
namespace _get_line {
namespace internal {
export {

// Recursively reads a line of arbitrary length.
std::string GetLine(std::ifstream& stream);

}  // export
}  // namespace internal

export using internal::GetLine;

}  // namespace _get_line
}  // namespace base
}  // namespace principia

#include "base/get_line_body.hpp"
