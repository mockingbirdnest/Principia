#pragma once

#include "ksp_plugin/part.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "ksp_plugin/frames.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_part {

using geometry::Velocity;
using physics::DegreesOfFreedom;

class MockPart : public Part {
 public:
  MockPart()
      : Part(/*part_id=*/666,
             /*mass=*/Mass(),
             /*degrees_of_freedom=*/
                 DegreesOfFreedom<Barycentric>(Barycentric::origin,
                                               Velocity<Barycentric>()),
             /*deletion_callback=*/nullptr) {}

  MOCK_CONST_METHOD0(part_id, PartId());
  MOCK_CONST_METHOD0(mass, Mass const&());
  MOCK_CONST_METHOD0(intrinsic_force, Vector<Force, Barycentric> const&());
  MOCK_CONST_METHOD0(degrees_of_freedom,
                     DegreesOfFreedom<Barycentric> const&());
  MOCK_CONST_METHOD0(tail, DiscreteTrajectory<Barycentric> const&());
  MOCK_CONST_METHOD0(tail_is_authoritative, bool());
};

}  // namespace internal_part
}  // namespace ksp_plugin
}  // namespace principia
