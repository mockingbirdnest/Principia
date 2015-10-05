#pragma once

#include "physics/body_centered_non_rotating_dynamic_frame.hpp"

#include "geometry/identity.hpp"

namespace principia {

using geometry::Identity;

namespace physics {

template<typename InertialFrame, typename ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
BodyCentredNonRotatingDynamicFrame(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    not_null<MassiveBody const*> const centre) 
    : ephemeris_(ephemeris),
      centre_(centre),
      centre_trajectory_(ephemeris_->trajectory(centre_)) {}

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
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             AngularVelocity<InertialFrame>(),
             centre_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<ThisFrame, InertialFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
FromThisFrameAtTime(Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const centre_degrees_of_freedom =
      centre_trajectory_->EvaluateDegreesOfFreedom(t, &hint_);
  RigidTransformation<ThisFrame, InertialFrame> const
      rigid_transformation(ThisFrame::origin,
                           centre_degrees_of_freedom.position(),
                           Identity<ThisFrame, InertialFrame>());
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             AngularVelocity<InertialFrame>(),
             centre_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
GeometricAcceleration(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  Vector<Acceleration, InertialFrame> const acceleration_of_centre =
      ephemeris_->ComputeGravitationalAcceleration(centre_, t);
  Vector<Acceleration, InertialFrame> const acceleration_at_point =
      ephemeris_->ComputeGravitationalAcceleration(
          degrees_of_freedom.position(), t);
  return acceleration_at_point - acceleration_of_centre;
}

}  // namespace physics
}  // namespace principia
