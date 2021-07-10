#pragma once

#include <type_traits>

namespace principia {
namespace numerics {

// A constexpr implementation of the IEEE 754:2008 scaleB function.
template<typename SourceFormat,
         typename LogBFormat,
         typename = std::enable_if_t<std::is_floating_point_v<SourceFormat> &&
                                     std::is_integral_v<LogBFormat>>>
constexpr SourceFormat ScaleB(SourceFormat x, LogBFormat N);

}  // namespace numerics
}  // namespace principia

#include "numerics/scale_b_body.hpp"
