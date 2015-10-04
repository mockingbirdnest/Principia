#pragma once

#include "physics/body_centered_non_rotating_dynamic_frame.hpp"

namespace principia {
namespace physics {

template<typename InertialFrame, typename ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
BodyCentredNonRotatingDynamicFrame(
    ContinuousTrajectory<InertialFrame> const& centre_trajectory) {}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  return RigidMotion<InertialFrame, ThisFrame>();
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<ThisFrame, InertialFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
FromThisFrameAtTime(Instant const& t) const {
  return RigidMotion<ThisFrame, InertialFrame>();
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
GeometricAcceleration(Instant const& t,
                      DegreesOfFreedom<ThisFrame> const & degrees_of_freedom) const {
  return Vector<Acceleration, ThisFrame>();
}

}  // namespace physics
}  // namespace principia
