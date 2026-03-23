#pragma once

#include <concepts>

namespace principia {
namespace numerics {
namespace _scale_b {
namespace internal {

// A constexpr implementation of the IEEE 754:2008 scaleB function.
template<std::floating_point SourceFormat, std::integral LogBFormat>
constexpr SourceFormat ScaleB(SourceFormat x, LogBFormat N);

}  // namespace internal

using internal::ScaleB;

}  // namespace _scale_b
}  // namespace numerics
}  // namespace principia

#include "numerics/scale_b_body.hpp"
