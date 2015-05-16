#pragma once

namespace principia {
namespace quantities {

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
class QuantityOrDoubleSerializer {};

}  // namespace quantities
}  // namespace principia

#include "quantities/serialization_body.hpp"
