#pragma once

#include "base/get_line.hpp"

#include <string>

namespace principia {
namespace base {

namespace {
int constexpr kBufferSize = 100;
}  // namespace

// This implementation is in O(N^2) because of the fixed-size buffer.  We could
// double the buffer size at each recursive call to make it O(N).
std::string GetLine(not_null<std::ifstream*> const stream) {
  char buffer[kBufferSize];
  if (!stream->getline(&buffer[0], kBufferSize).eof() && stream->fail()) {
    stream->clear();
    return buffer + GetLine(stream);
  }
  return buffer;
}

}  // namespace base
}  // namespace principia
