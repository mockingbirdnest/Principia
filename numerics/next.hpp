#pragma once

#include <type_traits>

namespace principia {
namespace numerics {
namespace _next {
namespace internal {

// A constexpr implementation of the IEEE 754:2008 nextUp and nextDown
// functions.
template<typename SourceFormat,
         typename = std::enable_if_t<std::is_floating_point_v<SourceFormat>>>
constexpr SourceFormat NextUp(SourceFormat x);
template<typename SourceFormat,
         typename = std::enable_if_t<std::is_floating_point_v<SourceFormat>>>
constexpr SourceFormat NextDown(SourceFormat x);

}  // namespace internal

using internal::NextDown;
using internal::NextUp;

}  // namespace _next
}  // namespace numerics
}  // namespace principia

#include "numerics/next_body.hpp"
