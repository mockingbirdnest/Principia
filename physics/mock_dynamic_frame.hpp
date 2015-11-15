#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/rigid_motion.hpp"

namespace principia {

using geometry::Instant;
using geometry::Vector;
using quantities::Acceleration;

namespace physics {

template<typename InertialFrame, typename ThisFrame>
class MockDynamicFrame : public DynamicFrame<InertialFrame, ThisFrame> {
 public:
  MockDynamicFrame() {}

  MOCK_CONST_METHOD1_T(ToThisFrameAtTime,
                       RigidMotion<InertialFrame, ThisFrame>(Instant const& t));
  MOCK_CONST_METHOD1_T(FromThisFrameAtTime,
                       RigidMotion<ThisFrame, InertialFrame>(Instant const& t));
  MOCK_CONST_METHOD2_T(
      GeometricAcceleration,
      Vector<Acceleration, ThisFrame>(
          Instant const& t,
          DegreesOfFreedom<ThisFrame> const& degrees_of_freedom));

  using Rot = Rotation<Frenet<ThisFrame>, ThisFrame>;

  MOCK_CONST_METHOD2_T(
      FrenetFrame,
      Rot(Instant const& t,
          DegreesOfFreedom<ThisFrame> const& degrees_of_freedom));
};

}  // namespace physics
}  // namespace principia

