#pragma once

#include "base/get_line.hpp"

#include <memory>
#include <string>

namespace principia {
namespace base {
namespace _get_line {
namespace internal {

constexpr int buffer_size = 200;

std::string GetLineWithSize(std::size_t const size, std::ifstream& stream) {
  std::unique_ptr<char[]> buffer(new char[size]);
  if (!stream.getline(&buffer[0], size).eof() && stream.fail()) {
    stream.clear();
    std::string string_buffer(buffer.get());
    string_buffer += GetLineWithSize(2 * size, stream);
    return string_buffer;
  } else {
    std::string string_buffer(buffer.get());
    return string_buffer;
  }
}

std::string GetLine(std::ifstream& stream) {
  return GetLineWithSize(buffer_size, stream);
}

}  // namespace internal
}  // namespace _get_line
}  // namespace base
}  // namespace principia
