#pragma once

#include <string>

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {

template<typename T>
T ParseQuantity(std::string const& s);

template<typename T>
T ParseUnit(std::string const& s);

}  // namespace quantities
}  // namespace principia

#include "quantities/parser_body.hpp"
