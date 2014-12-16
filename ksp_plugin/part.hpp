#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::Position;
using principia::geometry::Vector;
using principia::geometry::Velocity;
using principia::physics::DegreesOfFreedom;
using principia::quantities::Acceleration;
using principia::quantities::Mass;

namespace principia {
namespace ksp_plugin {

// Corresponds to KSP's |Part.flightID|, *not* to |Part.uid|.  C#'s |uint|
// corresponds to |uint32_t|.
using PartId = uint32_t;

// Represents a KSP part.
template<typename Frame>
struct Part {
  Part(DegreesOfFreedom<Frame> const& degrees_of_freedom,
       Mass const& mass,
       Vector<Acceleration, Frame> const&
           gravitational_acceleration_to_be_applied_by_ksp);

  DegreesOfFreedom<Frame> degrees_of_freedom;
  Mass mass;
  Vector<Acceleration, Frame>
      gravitational_acceleration_to_be_applied_by_ksp;
  // TODO(egg): we may want to keep track of the moment of inertia, angular
  // momentum, etc.
};

template<typename Frame>
std::ostream& operator<<(std::ostream& out, Part<Frame> const& part);

}  // namespace ksp_plugin
}  // namespace principia

#include "part_body.hpp"
