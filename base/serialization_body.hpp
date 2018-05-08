#pragma once

#include "base/serialization.hpp"

namespace principia {
namespace base {

inline UniqueBytes SerializeAsBytes(
    google::protobuf::MessageLite const& message) {
  UniqueBytes bytes(message.ByteSizeLong());
  message.SerializeToArray(bytes.data.get(), bytes.size);
  return std::move(bytes);
}

}  // namespace base
}  // namespace principia
