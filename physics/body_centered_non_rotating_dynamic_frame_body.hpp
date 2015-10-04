#pragma once

#include "physics/body_centered_non_rotating_dynamic_frame.hpp"

#include "geometry/identity.hpp"

namespace principia {

using geometry::Identity;

namespace physics {

template<typename InertialFrame, typename ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
BodyCentredNonRotatingDynamicFrame(
    not_null<ContinuousTrajectory<InertialFrame> const*> const
        centre_trajectory) 
    : centre_trajectory_(centre_trajectory) {}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const centre_degrees_of_freedom =
      centre_trajectory_->EvaluateDegreesOfFreedom(t, &hint_);
  RigidTransformation<InertialFrame, ThisFrame> const
      rigid_transformation(centre_degrees_of_freedom.position(),
                           ThisFrame::origin,
                           Identity<InertialFrame, ThisFrame>());
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
