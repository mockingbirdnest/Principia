#pragma once

namespace principia {
namespace geometry {

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
class DoubleOrQuantityOrMultivectorSerializer {};

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
class PointOrMultivectorSerializer {};

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
class QuantityOrMultivectorSerializer {};

}  // namespace geometry
}  // namespace principia

#include "geometry/serialization_body.hpp"
