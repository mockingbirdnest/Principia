export module principia.base.serialization;

import "google/protobuf/message_lite.h";

import principia.base.array;
import principia.base.serialization;

namespace principia {
namespace base {
namespace _serialization {
namespace internal {

using namespace principia::base::_array;

export {

inline UniqueArray<std::uint8_t> SerializeAsBytes(
    google::protobuf::MessageLite const& message);

template<typename Message>
Message ParseFromBytes(Array<std::uint8_t const> bytes);

}  // export
}  // namespace internal

export using internal::ParseFromBytes;
export using internal::SerializeAsBytes;

}  // namespace _serialization
}  // namespace base
}  // namespace principia

#include "base/serialization_body.hpp"
