#pragma once

#include <type_traits>

namespace principia {
namespace numerics {
namespace _log_b {
namespace internal {

// A constexpr implementation of the IEEE 754:2008 logB function.
// Uses sourceFormat as logBFormat, which makes it easy to cleanly handle NaN,
// infinity, and 0.
template<typename SourceFormat,
         typename = std::enable_if_t<std::is_floating_point_v<SourceFormat>>>
constexpr SourceFormat LogB(SourceFormat x);

}  // namespace internal

using internal::LogB;

}  // namespace _log_b
}  // namespace numerics
}  // namespace principia

namespace principia::numerics {
using namespace principia::numerics::_log_b;
}  // namespace principia::numerics

#include "numerics/log_b_body.hpp"
