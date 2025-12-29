#pragma once

#include <utility>
#include <vector>

#include "nanobenchmarks/nanobenchmark.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _microarchitectures {
namespace internal {

using namespace principia::nanobenchmarks::_nanobenchmark;

template<typename Value, typename Argument>
using NanobenchmarkAndCycles =
    std::pair<Nanobenchmark<Value, Argument> const*, int>;

// The returned vector is sorted by nanobenchmark name.
template<typename Value, typename Argument>
std::vector<NanobenchmarkAndCycles<Value, Argument>> const&
ReferenceCycleCounts();

}  // namespace internal

using internal::ReferenceCycleCounts;

}  // namespace _microarchitectures
}  // namespace nanobenchmarks
}  // namespace principia

#include "nanobenchmarks/microarchitectures_body.hpp"
