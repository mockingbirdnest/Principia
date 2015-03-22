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

inline Bytes Bytes::New(std::int64_t const size) {
  return Bytes(new std::uint8_t[size], size);
}

inline Bytes Bytes::New(std::string s) {
  Bytes result = New(s.size());
  std::memcpy(result.data, s.c_str(), s.size());
}

inline void Bytes::Delete(Bytes const bytes) {
  delete[] bytes.data;
}

inline bool operator==(Bytes left, Bytes right) {
  if (left.size != right.size) {
    return false;
  }
  return std::memcmp(left.data, right.data, right.size) == 0;
}

inline void BytesDeleter::operator()(Bytes* const bytes) const {
  Bytes::Delete(*bytes);
  delete bytes;
}

}  // namespace base
}  // namespace principia
