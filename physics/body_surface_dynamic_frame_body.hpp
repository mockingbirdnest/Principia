
#pragma once

#include "physics/body_surface_dynamic_frame.hpp"

#include "geometry/identity.hpp"

namespace principia {
namespace physics {
namespace internal_body_surface_dynamic_frame {

using geometry::AngularVelocity;
using geometry::Identity;

template<typename InertialFrame, typename ThisFrame>
BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::
BodySurfaceDynamicFrame(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    not_null<RotatingBody<InertialFrame> const*> const centre)
    : ephemeris_(ephemeris),
      centre_(centre),
      centre_trajectory_(ephemeris_->trajectory(centre_)) {}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const centre_degrees_of_freedom =
      centre_trajectory_->EvaluateDegreesOfFreedom(t, &hint_);
  //TODO(phl):Change
  RigidTransformation<InertialFrame, ThisFrame> const
      rigid_transformation(centre_degrees_of_freedom.position(),
                           ThisFrame::origin,
                           Identity<InertialFrame, ThisFrame>().Forget());
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             AngularVelocity<InertialFrame>(),
             centre_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
void BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::
WriteToMessage(not_null<serialization::DynamicFrame*> const message) const {
  message->MutableExtension(
      serialization::BodySurfaceDynamicFrame::extension)->set_centre(
          ephemeris_->serialization_index_for_body(centre_));
}

template<typename InertialFrame, typename ThisFrame>
not_null<std::unique_ptr<
    BodySurfaceDynamicFrame<InertialFrame, ThisFrame>>>
BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::BodySurfaceDynamicFrame const& message) {
  return std::make_unique<BodySurfaceDynamicFrame>(
             ephemeris,
             ephemeris->body_for_serialization_index(message.centre()));
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, InertialFrame>
BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::
    GravitationalAcceleration(Instant const& t,
                              Position<InertialFrame> const& q) const {
  return ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(q, t);
}

template<typename InertialFrame, typename ThisFrame>
AcceleratedRigidMotion<InertialFrame, ThisFrame>
BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::MotionOfThisFrame(
    Instant const& t) const {
  //TODO(phl):Change
  return AcceleratedRigidMotion<InertialFrame, ThisFrame>(
             ToThisFrameAtTime(t),
             /*angular_acceleration_of_to_frame=*/{},
             /*acceleration_of_to_frame_origin=*/ephemeris_->
                 ComputeGravitationalAccelerationOnMassiveBody(centre_, t));
}

}  // namespace internal_body_surface_dynamic_frame
}  // namespace physics
}  // namespace principia
