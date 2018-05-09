#pragma once

#include "base/serialization.hpp"

namespace principia {
namespace base {

inline UniqueArray<std::uint8_t> SerializeAsBytes(
    google::protobuf::MessageLite const& message) {
  UniqueArray<std::uint8_t> bytes(message.ByteSizeLong());
  message.SerializeToArray(bytes.data.get(), bytes.size);
  return std::move(bytes);
}

}  // namespace base
}  // namespace principia
