#pragma once

#include "geometry/frame.hpp"
#include "geometry/rotation.hpp"
#include "physics/rigid_motion.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using geometry::Rotation;

namespace physics {

// The Frenet frame of a free fall trajectory in |Frame|.
// TODO(egg): this should actually depend on its template parameter somehow.
template<typename Frame>
using Frenet = geometry::Frame<serialization::Frame::PhysicsTag,
                                    serialization::Frame::FRENET,
                                    false /*frame_is_inertial*/>;

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

  // The acceleration due to the non-inertial motion of |ThisFrame| and gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is |GeometricAcceleration|.
  virtual Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const = 0;

  // The definition of the Frenet frame of a free fall trajectory in |ThisFrame|
  // with the given |degrees_of_freedom| at instant |t|.
  virtual Rotation<Frenet<ThisFrame>, ThisFrame> FrenetFrame(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const;
};

}  // namespace physics
}  // namespace principia

#include "physics/dynamic_frame_body.hpp"
