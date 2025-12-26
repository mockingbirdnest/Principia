#pragma once

#include "geometry/instant.hpp"
#include "geometry/space.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _dependencies {
namespace internal {

using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;

template<typename Value_, typename Argument_>
struct Dependencies;

template<>
struct Dependencies<double, double> {
  static double ProduceArgument(double x);
  static double ConsumeValue(double value);
};

template<typename Frame>
struct Dependencies<Displacement<Frame>, Instant> {
  static Instant ProduceArgument(double x);
  static double ConsumeValue(Displacement<Frame> const& value);
};

}  // namespace internal
}  // namespace _dependencies
}  // namespace nanobenchmarks
}  // namespace principia

#include "nanobenchmarks/dependencies_body.hpp"