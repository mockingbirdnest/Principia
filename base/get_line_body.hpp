
#pragma once

#include "base/get_line.hpp"

#include <string>

namespace principia {
namespace base {
namespace internal_get_line {

constexpr int buffer_size = 200;

std::string GetLineWithSize(std::size_t const size,
                            not_null<std::ifstream*> const stream) {
  std::unique_ptr<char[]> buffer(new char[size]);
  if (!stream->getline(&buffer[0], size).eof() && stream->fail()) {
    stream->clear();
    std::string string_buffer(buffer.get());
    string_buffer += GetLineWithSize(2 * size, stream);
    return std::move(string_buffer);
  } else {
    std::string string_buffer(buffer.get());
    return std::move(string_buffer);
  }
}

std::string GetLine(not_null<std::ifstream*> const stream) {
  return GetLineWithSize(buffer_size, stream);
}

}  // namespace internal_get_line
}  // namespace base
}  // namespace principia
