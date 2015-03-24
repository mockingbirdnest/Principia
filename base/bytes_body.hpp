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

UniqueBytes::UniqueBytes() : size(0) {}

UniqueBytes::UniqueBytes(std::int64_t const size)
    : data(std::make_unique<std::uint8_t[]>(static_cast<size_t>(size))),
      size(size) {}

UniqueBytes::UniqueBytes(std::unique_ptr<std::uint8_t[]> data,
                         std::int64_t const size)
    : data(std::move(data)),
      size(size) {}

inline bool operator==(Bytes left, Bytes right) {
  if (left.size != right.size) {
    return false;
  }
  return std::memcmp(left.data,
                     right.data,
                     static_cast<size_t>(right.size)) == 0;
}

inline bool operator==(UniqueBytes left, UniqueBytes right) {
  if (left.size != right.size) {
    return false;
  }
  return std::memcmp(left.data.get(),
                     right.data.get(),
                     static_cast<size_t>(right.size)) == 0;
}

}  // namespace base
}  // namespace principia
