#pragma once

#include <map>
#include <memory>

#include "ksp_plugin/frames.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

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
class Part {
 public:
  Part(DegreesOfFreedom<Frame> const& degrees_of_freedom,
       Mass const& mass,
       Vector<Acceleration, Frame> const&
           gravitational_acceleration_to_be_applied_by_ksp);

  //Part(Part&&) = delete;
  //Part& operator=(Part&&) = delete;
  Part(Part const&) = delete;
  Part& operator=(Part const&) = delete;

  DegreesOfFreedom<Frame> const& degrees_of_freedom() const;
  Mass const& mass() const;
  Vector<Acceleration, Frame> const&
      gravitational_acceleration_to_be_applied_by_ksp() const;

  void WriteToMessage(not_null<serialization::Part*> const message) const;
  static Part ReadFromMessage(serialization::Part const& message);

 private:
  DegreesOfFreedom<Frame> degrees_of_freedom_;
  Mass mass_;
  Vector<Acceleration, Frame>
      gravitational_acceleration_to_be_applied_by_ksp_;
  // TODO(egg): we may want to keep track of the moment of inertia, angular
  // momentum, etc.
};

template<typename Frame>
std::ostream& operator<<(std::ostream& out, Part<Frame> const& part);

using PartIdToOwnedPart = std::map<PartId,
                                   not_null<std::unique_ptr<Part<World>>>>;
using IdAndOwnedPart = PartIdToOwnedPart::value_type;

}  // namespace ksp_plugin
}  // namespace principia

#include "part_body.hpp"
