#pragma once

#include "base/get_line.hpp"

#include <string>

namespace principia {
namespace base {

namespace {

int constexpr kBufferSize = 100;

std::string GetLineWithSize(std::size_t const size,
                            not_null<std::ifstream*> const stream) {
  char* buffer = new char[size];
  if (!stream->getline(buffer, size).eof() && stream->fail()) {
    stream->clear();
    std::string string_buffer(buffer);
    string_buffer += GetLineWithSize(2 * size, stream);
    return std::move(string_buffer);
  }
  return buffer;
}

}  // namespace

std::string GetLine(not_null<std::ifstream*> const stream) {
  return GetLineWithSize(kBufferSize, stream);
}

}  // namespace base
}  // namespace principia
