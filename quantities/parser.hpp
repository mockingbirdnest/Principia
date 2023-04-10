#pragma once

#include <string>

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace _parser {
namespace internal {

// This function parses the following grammar:
//   quantity            ⩴ double quotient_unit
//   quotient_unit       ⩴ quotient_unit / exponentiation_unit
//                       | / exponentiation_unit
//                       | product_unit
//   product_unit        ⩴ exponentiation_unit [product_unit]
//   exponentiation_unit ⩴ unit [^ exponent]
//   exponent            ⩴ signed_integer
template<typename Q>
Q ParseQuantity(std::string const& s);

}  // namespace internal

using internal::ParseQuantity;

}  // namespace _parser
}  // namespace quantities
}  // namespace principia

#include "quantities/parser_body.hpp"
