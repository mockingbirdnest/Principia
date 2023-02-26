#pragma once

#include "base/not_constructible.hpp"

namespace principia {
namespace quantities {
namespace internal_serialization {

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

}  // namespace internal_serialization

using internal_serialization::DoubleOrQuantitySerializer;

}  // namespace quantities
}  // namespace principia

#include "quantities/serialization_body.hpp"
