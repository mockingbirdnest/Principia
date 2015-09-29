#pragma once

#include "physics/rigid_motion.hpp"

namespace principia {
namespace physics {

// The definition of a reference frame |ThisFrame| in arbitrary motion with
// respect to the inertial reference frame |InertialFrame|.
template<typename InertialFrame, typename ThisFrame>
class DynamicFrame {
  static_assert(InertialFrame::is_inertial, "InertialFrame must be inertial");

 public:
  virtual ~DynamicFrame() = default;
  virtual RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const = 0;
  virtual RigidMotion<ThisFrame, InertialFrame> FromThisFrameAtTime(
      Instant const& t) const = 0;

  // The acceleration due to the non-inertial motion of |Frame| and gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is |GeometricAcceleration|.
  virtual Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const = 0;
};

}  // namespace physics
}  // namespace principia
