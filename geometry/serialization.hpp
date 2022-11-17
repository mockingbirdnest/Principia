#pragma once

#include "base/not_constructible.hpp"

namespace principia {
namespace geometry {
namespace internal_serialization {

using base::not_constructible;

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

}  // namespace internal_serialization

using internal_serialization::DoubleOrQuantityOrPointOrMultivectorSerializer;
using internal_serialization::DoubleOrQuantityOrMultivectorSerializer;
using internal_serialization::PointOrMultivectorSerializer;
using internal_serialization::QuantityOrMultivectorSerializer;

}  // namespace geometry
}  // namespace principia

#include "geometry/serialization_body.hpp"
