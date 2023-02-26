#pragma once

#include "base/array.hpp"

#include "google/protobuf/message_lite.h"

namespace principia {
namespace base {
namespace _serialization {
namespace internal {

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

namespace principia::base {
using namespace principia::base::_serialization;
}  // namespace principia::base

#include "base/serialization_body.hpp"
