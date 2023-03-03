#pragma once

#include "base/not_constructible.hpp"

namespace principia {
namespace quantities {
namespace _serialization {
namespace internal {

using namespace principia::base::_not_constructible;

// A helper class that serializes a |Quantity| or a |double| to a protobuf
// structure like:
//
// message Message {
//   oneof message {
//     double double = 1;
//     Quantity quantity = 2;
//   }
// }
template<typename T, typename Message>
struct DoubleOrQuantitySerializer : not_constructible {};

}  // namespace internal

using internal::DoubleOrQuantitySerializer;

}  // namespace _serialization
}  // namespace quantities
}  // namespace principia

namespace principia::quantities {
using namespace principia::quantities::_serialization;
}  // namespace principia::quantities

#include "quantities/serialization_body.hpp"
