#pragma once

#include "base/not_constructible.hpp"

namespace principia {
namespace geometry {
namespace _serialization {
namespace internal {

using namespace principia::base::_not_constructible;

// A helper class that serializes a |double|, a |Quantity|, a |Point| or a
// |Multivector| to a protobuf structure like:
//
// message Message {
//   oneof message {
//     double double = 1;
//     Quantity quantity = 2;
//     Point point = 3;
//     Multivector multivector = 4;
//   }
// }
template<typename T, typename Message>
struct DoubleOrQuantityOrPointOrMultivectorSerializer : not_constructible {};

// A helper class that serializes a |double|, a |Quantity| or a |Multivector|
// to a protobuf structure like:
//
// message Message {
//   oneof message {
//     double double = 1;
//     Quantity quantity = 2;
//     Multivector multivector = 3;
//   }
// }
template<typename T, typename Message>
struct DoubleOrQuantityOrMultivectorSerializer : not_constructible {};

// A helper class that serializes a |Point| or a |Multivector| to a protobuf
// structure like:
//
// message Message {
//   oneof message {
//     Point point = 1;
//     Multivector multivector = 2;
//   }
// }
template<typename T, typename Message>
struct PointOrMultivectorSerializer : not_constructible {};

// A helper class that serializes a |Quantity| or a |Multivector| to a protobuf
// structure like:
//
// message Message {
//   oneof message {
//     Quantity quantity = 1;
//     Multivector multivector = 2;
//   }
// }
template<typename T, typename Message>
struct QuantityOrMultivectorSerializer : not_constructible {};

}  // namespace internal

using internal::DoubleOrQuantityOrPointOrMultivectorSerializer;
using internal::DoubleOrQuantityOrMultivectorSerializer;
using internal::PointOrMultivectorSerializer;
using internal::QuantityOrMultivectorSerializer;

}  // namespace _serialization
}  // namespace geometry
}  // namespace principia

#include "geometry/serialization_body.hpp"
