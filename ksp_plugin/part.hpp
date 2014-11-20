#pragma once

#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::Position;
using principia::geometry::Velocity;
using principia::quantities::Mass;

namespace principia {
namespace ksp_plugin {

// Represents a KSP part.
template<typename Frame>
struct Part {
  Position<Frame> position;
  Velocity<Frame> velocity;
  Mass mass;
  // TODO(egg): we may want to keep track of the moment of inertia, angular
  // momentum, etc.
};

}  // namespace ksp_plugin
}  // namespace principia
