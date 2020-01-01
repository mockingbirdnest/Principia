
#pragma once

#include "physics/body_surface_dynamic_frame.hpp"

#include "base/not_null.hpp"
#include "geometry/rotation.hpp"

namespace principia {
namespace physics {
namespace internal_body_surface_dynamic_frame {

using base::check_not_null;
using base::dynamic_cast_not_null;
using geometry::AngularVelocity;
using geometry::OrthogonalMap;
using geometry::Rotation;
using quantities::Variation;

template<typename InertialFrame, typename ThisFrame>
BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::
BodySurfaceDynamicFrame(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    not_null<RotatingBody<InertialFrame> const*> const centre)
    : ephemeris_(ephemeris),
      centre_(centre),
      centre_trajectory_(ephemeris_->trajectory(centre_)) {}

template<typename InertialFrame, typename ThisFrame>
not_null<RotatingBody<InertialFrame> const*>
BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::centre() const {
  return centre_;
}

template<typename InertialFrame, typename ThisFrame>
Instant BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::t_min()
    const {
  return centre_trajectory_->t_min();
}

template<typename InertialFrame, typename ThisFrame>
Instant BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::t_max()
    const {
  return centre_trajectory_->t_max();
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BodySurfaceDynamicFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const centre_degrees_of_freedom =
      centre_trajectory_->EvaluateDegreesOfFreedom(t);

  Rotation<InertialFrame, ThisFrame> rotation =
      centre_->template ToSurfaceFrame<ThisFrame>(t);
  AngularVelocity<InertialFrame> angular_velocity = centre_->angular_velocity();
  RigidTransformation<InertialFrame, ThisFrame> const
      rigid_transformation(centre_degrees_of_freedom.position(),
                           ThisFrame::origin,
                           rotation.Forget<OrthogonalMap>());
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             angular_velocity,
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
             dynamic_cast_not_null<RotatingBody<InertialFrame> const*>(
                 ephemeris->body_for_serialization_index(message.centre())));
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
  DegreesOfFreedom<InertialFrame> const centre_degrees_of_freedom =
      centre_trajectory_->EvaluateDegreesOfFreedom(t);
  Vector<Acceleration, InertialFrame> const centre_acceleration =
      ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(centre_, t);

  auto const to_this_frame = ToThisFrameAtTime(t);

  Variation<AngularVelocity<InertialFrame>> const
      angular_acceleration_of_to_frame;
  Vector<Acceleration, InertialFrame> const& acceleration_of_to_frame_origin =
      centre_acceleration;
  return AcceleratedRigidMotion<InertialFrame, ThisFrame>(
             to_this_frame,
             angular_acceleration_of_to_frame,
             acceleration_of_to_frame_origin);
}

}  // namespace internal_body_surface_dynamic_frame
}  // namespace physics
}  // namespace principia
