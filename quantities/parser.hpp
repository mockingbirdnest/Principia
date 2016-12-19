
#pragma once

#include <string>

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace internal_parser {

template<typename T>
T ParseQuantity(std::string const& s);

template<typename T>
T ParseUnit(std::string const& s);

}  // namespace internal_parser

using internal_parser::ParseQuantity;
using internal_parser::ParseUnit;

}  // namespace quantities
}  // namespace principia

#include "quantities/parser_body.hpp"
