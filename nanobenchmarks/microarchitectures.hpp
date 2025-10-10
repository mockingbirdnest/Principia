#pragma once

#include <map>

#include "nanobenchmarks/function_registry.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _microarchitectures {
namespace internal {

using namespace principia::nanobenchmarks::_function_registry;

std::vector<std::pair<Nanobenchmark const*, int>> const& ReferenceCycleCounts();

}  // namespace internal

using internal::ReferenceCycleCounts;

}  // namespace _microarchitectures
}  // namespace nanobenchmarks
}  // namespace principia
