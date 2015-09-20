#pragma once

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {

template<typename Q>
Q ParseQuantity(std::string const& s);

template<typename Q>
Q ParseUnit(std::string const& s);

}  // namespace quantities
}  // namespace principia

#include "quantities/parser_body.hpp"
