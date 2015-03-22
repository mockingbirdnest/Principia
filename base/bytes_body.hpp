#pragma once

#include "base/bytes.hpp"

#include "glog/logging.h"

namespace principia {
namespace base {

inline Bytes::Bytes()
    : data(nullptr), size(0) {}

inline Bytes::Bytes(std::uint8_t* const data,
                    std::int64_t const size)
    : data(data), size(size) {}

inline void Bytes::CheckNotNull() const {
  CHECK_NOTNULL(data);
}

Bytes Bytes::New(std::int64_t const size) {
  return Bytes(new std::uint8_t[size], size);
}

Bytes Bytes::New(std::string s) {
  Bytes result = New(s.size());
  std::memcpy(result.data, s.c_str(), s.size());
}

void Bytes::Delete(Bytes const bytes) {
  delete[] bytes.data;
}

}  // namespace base
}  // namespace principia
