
#pragma once

#include <string>

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace internal_parser {

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

}  // namespace internal_parser

using internal_parser::ParseQuantity;

}  // namespace quantities
}  // namespace principia

#include "quantities/parser_body.hpp"
