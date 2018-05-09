#pragma once

#include "base/array.hpp"

#include "google/protobuf/message_lite.h"

namespace principia {
namespace base {

UniqueArray<std::uint8_t> SerializeAsBytes(
    google::protobuf::MessageLite const& message);

template<typename Message>
Message ParseFromBytes(Array<std::uint8_t const> bytes);

}  // namespace base
}  // namespace principia

#include "base/serialization_body.hpp"
