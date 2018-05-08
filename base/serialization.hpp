#pragma once

#include "base/array.hpp"

#include "google/protobuf/message_lite.h"

namespace principia {
namespace base {

inline UniqueBytes SerializeAsBytes(
    google::protobuf::MessageLite const& message);

}  // namespace base
}  // namespace principia

#include "base/serialization_body.hpp"
