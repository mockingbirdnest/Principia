#pragma once

#include "base/array.hpp"

#include "google/protobuf/message_lite.h"

namespace principia {
namespace base {

inline UniqueArray<std::uint8_t> SerializeAsBytes(
    google::protobuf::MessageLite const& message);

}  // namespace base
}  // namespace principia

#include "base/serialization_body.hpp"
