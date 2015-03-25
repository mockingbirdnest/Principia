#pragma once

#include "base/bytes.hpp"

#include "glog/logging.h"

namespace principia {
namespace base {

inline Bytes::Bytes()
    : data(nullptr), size(0) {}

template<typename T>
Bytes::Bytes(std::uint8_t* const data,
             T const size)
    : data(data), size(static_cast<std::int64_t>(size)) {}

inline UniqueBytes::UniqueBytes() : size(0) {}

template<typename T>
UniqueBytes::UniqueBytes(T const size)
    : data(new std::uint8_t[static_cast<size_t>(size)]),
      size(static_cast<std::int64_t>(size)) {}

template<typename T>
UniqueBytes::UniqueBytes(std::unique_ptr<std::uint8_t[]> data,
                         T const size)
    : data(data.release()),
      size(static_cast<std::int64_t>(size)) {}

inline UniqueBytes::~UniqueBytes() {
  delete[] data;
  data = nullptr;
}

inline bool operator==(Bytes left, Bytes right) {
  if (left.size != right.size) {
    return false;
  }
  return std::memcmp(left.data,
                     right.data,
                     static_cast<size_t>(right.size)) == 0;
}

inline bool operator==(Bytes left, UniqueBytes right) {
  return left == Bytes(right.data, right.size);
}

inline bool operator==(UniqueBytes left, Bytes right) {
  return Bytes(left.data, left.size) == right;
}

inline bool operator==(UniqueBytes left, UniqueBytes right) {
  return Bytes(left.data, left.size) == Bytes(right.data, right.size);
}

}  // namespace base
}  // namespace principia
