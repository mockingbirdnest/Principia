#pragma once

#include <algorithm>

#include "base/hexadecimal.hpp"
#include "base/serialization.hpp"

namespace principia {
namespace base {

UniqueArray<std::uint8_t> SerializeAsBytes(
    google::protobuf::MessageLite const& message) {
  UniqueArray<std::uint8_t> bytes(message.ByteSizeLong());
  message.SerializeToArray(bytes.data.get(), bytes.size);
  return std::move(bytes);
}

template<typename Message>
Message ParseFromBytes(Array<std::uint8_t const> bytes) {
  Message message;
  CHECK(message.ParseFromArray(bytes.data, bytes.size));
  return message;
}

}  // namespace base
}  // namespace principia
