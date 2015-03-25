#pragma once

#include "base/bytes.hpp"

#include "glog/logging.h"

namespace principia {
namespace base {

inline Bytes::Bytes()
    : data(nullptr), size(0) {}

inline Bytes::Bytes(UniqueBytes const& bytes)
    : data(bytes.data.get()), size(bytes.size) {}

template<typename T, typename>
Bytes::Bytes(std::uint8_t* const data, T const size)
    : data(data), size(static_cast<std::int64_t>(size)) {}

inline UniqueBytes::UniqueBytes() : size(0) {}

template<typename T, typename>
UniqueBytes::UniqueBytes(T const size)
    : data(new std::uint8_t[static_cast<size_t>(size)]),
      size(static_cast<std::int64_t>(size)) {}

template<typename T, typename>
UniqueBytes::UniqueBytes(std::unique_ptr<std::uint8_t[]> data, T const size)
    : data(data.release()),
      size(static_cast<std::int64_t>(size)) {}

inline UniqueBytes::UniqueBytes(UniqueBytes&& bytes)
    : data(std::move(bytes.data)), size(bytes.size) {}

inline UniqueBytes& UniqueBytes::operator=(UniqueBytes&& bytes) {
  data = std::move(bytes.data);
  size = bytes.size;
  return *this;
}

inline bool operator==(Bytes left, Bytes right) {
  if (left.size != right.size) {
    return false;
  }
  return std::memcmp(left.data,
                     right.data,
                     static_cast<size_t>(right.size)) == 0;
}

inline bool operator==(Bytes left, UniqueBytes const& right) {
  return left == Bytes(right);
}

inline bool operator==(UniqueBytes const& left, Bytes right) {
  return Bytes(left) == right;
}

inline bool operator==(UniqueBytes const& left, UniqueBytes const& right) {
  return Bytes(left) == Bytes(right);
}

}  // namespace base
}  // namespace principia
