#pragma once

#include <concepts>

namespace principia {
namespace numerics {
namespace _log_b {
namespace internal {

// A constexpr implementation of the IEEE 754:2008 logB function.
// Uses sourceFormat as logBFormat, which makes it easy to cleanly handle NaN,
// infinity, and 0.
template<std::floating_point SourceFormat>
constexpr SourceFormat LogB(SourceFormat x);

}  // namespace internal

using internal::LogB;

}  // namespace _log_b
}  // namespace numerics
}  // namespace principia

#include "numerics/log_b_body.hpp"
