#pragma once

#include "base/array.hpp"

#include "google/protobuf/message_lite.h"

namespace principia {
namespace base {
namespace _serialization {
namespace internal {

using namespace principia::base::_array;

inline UniqueArray<std::uint8_t> SerializeAsBytes(
    google::protobuf::MessageLite const& message);

template<typename Message>
Message ParseFromBytes(Array<std::uint8_t const> bytes);

}  // namespace internal

using internal::ParseFromBytes;
using internal::SerializeAsBytes;

}  // namespace _serialization
}  // namespace base
}  // namespace principia

#include "base/serialization_body.hpp"
