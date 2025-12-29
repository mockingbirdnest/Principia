#pragma once

#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _dependencies {
namespace internal {

using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;

// NOTE: Add a specialization of `Dependencies` for each pair of `Value,
// Argument` used in nanobenchmarks.

// A concrete frame useful for declaring actual nanobenchmarks.
using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

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
  static Instant ProduceArgument(double x);
  static double ConsumeValue(Displacement<Frame> const& value);

 private:
  constexpr static Instant t0_;
};

template<typename Frame>
class Dependencies<RelativeDegreesOfFreedom<Frame>, Instant> {
 public:
  static Instant ProduceArgument(double x);
  static double ConsumeValue(RelativeDegreesOfFreedom<Frame> const& value);

 private:
  constexpr static Instant t0_;
};

}  // namespace internal

using internal::Dependencies;
using internal::World;

}  // namespace _dependencies
}  // namespace nanobenchmarks
}  // namespace principia

#include "nanobenchmarks/dependencies_body.hpp"
