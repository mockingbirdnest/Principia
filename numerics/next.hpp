#pragma once

#include <concepts>

namespace principia {
namespace numerics {
namespace _next {
namespace internal {

// A constexpr implementation of the IEEE 754:2008 nextUp and nextDown
// functions.
template<std::floating_point SourceFormat>
constexpr SourceFormat NextUp(SourceFormat x);
template<std::floating_point SourceFormat>
constexpr SourceFormat NextDown(SourceFormat x);

}  // namespace internal

using internal::NextDown;
using internal::NextUp;

}  // namespace _next
}  // namespace numerics
}  // namespace principia

#include "numerics/next_body.hpp"
