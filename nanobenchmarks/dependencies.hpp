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
class Dependencies;

template<>
class Dependencies<double, double> {
 public:
  static double ProduceArgument(double x);
  static double ConsumeValue(double value);
};

template<typename Frame>
class Dependencies<Displacement<Frame>, Instant> {
 public:
  static constexpr int expected_cycles = 16;//TODO(phl):Fix
  static Instant ProduceArgument(double x);
  static Displacement<Frame> Run(Instant argument);
  static double ConsumeValue(Displacement<Frame> const& value);

 private:
  constexpr static Instant t0_;
};

}  // namespace internal
}  // namespace _dependencies
}  // namespace nanobenchmarks
}  // namespace principia

#include "nanobenchmarks/dependencies_body.hpp"
