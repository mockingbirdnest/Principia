#pragma once

namespace principia {
namespace geometry {

// A helper class that serializes a |double| or a |Quantity| to a protobuf
// structure like:
//
// message Message {
//   oneof message {
//     double double = 1;
//     Quantity quantity = 2;
//   }
// }
template<typename T, typename Message>
class DoubleOrQuantitySerializer {};

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

}  // namespace geometry
}  // namespace principia

#include "geometry/serialization_body.hpp"
