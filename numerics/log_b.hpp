#pragma once

#include <type_traits>

namespace principia {
namespace numerics {

// A constexpr implementation of the IEEE 754:2008 logB function.
// Uses sourceFormat as logBFormat, which makes it easy to cleanly handle NaN,
// infinity, and 0.
template<typename SourceFormat,
         typename = std::enable_if_t<std::is_floating_point_v<SourceFormat>>>
constexpr SourceFormat LogB(SourceFormat const x);

}  // namespace numerics
}  // namespace principia

#include "numerics/log_b_body.hpp"
