#pragma once

#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::Position;
using principia::geometry::Velocity;
using principia::quantities::Mass;

namespace principia {
namespace ksp_plugin {

// Corresponds to KSP's |Part.flightID|, *not* to |Part.uid|.  C#'s |uint|
// corresponds to |uint32_t|.
using PartID = uint32_t;

// Represents a KSP part.
template<typename Frame>
struct Part {
  Position<Frame> position;
  Velocity<Frame> velocity;
  Mass mass;
  Vector<Acceleration, Frame> expected_ksp_acceleration;
  PartID id;
  // TODO(egg): we may want to keep track of the moment of inertia, angular
  // momentum, etc.
};

}  // namespace ksp_plugin
}  // namespace principia
