#pragma once

#include "physics/body_surface_reference_frame.hpp"

#include <utility>

#include "base/not_null.hpp"
#include "geometry/rotation.hpp"

namespace principia {
namespace physics {
namespace _body_surface_reference_frame {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_rotation;
using namespace principia::quantities::_named_quantities;

template<typename InertialFrame, typename ThisFrame>
BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::
BodySurfaceReferenceFrame(
    not_null<Ephemeris<InertialFrame> const*> ephemeris,
    not_null<RotatingBody<InertialFrame> const*> centre)
    : ephemeris_(std::move(ephemeris)),
      centre_(std::move(centre)),
      centre_trajectory_(ephemeris_->trajectory(centre_)) {}

template<typename InertialFrame, typename ThisFrame>
not_null<RotatingBody<InertialFrame> const*>
BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::centre() const {
  return centre_;
}

template<typename InertialFrame, typename ThisFrame>
Instant BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::t_min()
    const {
  return centre_trajectory_->t_min();
}

template<typename InertialFrame, typename ThisFrame>
Instant BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::t_max()
    const {
  return centre_trajectory_->t_max();
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const centre_degrees_of_freedom =
      centre_trajectory_->EvaluateDegreesOfFreedom(t);

  Rotation<InertialFrame, ThisFrame> rotation =
      centre_->template ToSurfaceFrame<ThisFrame>(t);
  AngularVelocity<InertialFrame> angular_velocity = centre_->angular_velocity();
  RigidTransformation<InertialFrame, ThisFrame> const
      rigid_transformation(centre_degrees_of_freedom.position(),
                           ThisFrame::origin,
                           rotation.template Forget<OrthogonalMap>());
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             angular_velocity,
             centre_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
void BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::
WriteToMessage(not_null<serialization::RigidReferenceFrame*> const message) const {
  message->MutableExtension(
      serialization::BodySurfaceReferenceFrame::extension)->set_centre(
          ephemeris_->serialization_index_for_body(centre_));
}

template<typename InertialFrame, typename ThisFrame>
not_null<std::unique_ptr<
    BodySurfaceReferenceFrame<InertialFrame, ThisFrame>>>
BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::BodySurfaceReferenceFrame const& message) {
  return std::make_unique<BodySurfaceReferenceFrame>(
             ephemeris,
             dynamic_cast_not_null<RotatingBody<InertialFrame> const*>(
                 ephemeris->body_for_serialization_index(message.centre())));
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, InertialFrame>
BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::
GravitationalAcceleration(Instant const& t,
                          Position<InertialFrame> const& q) const {
  return ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(q, t);
}

template<typename InertialFrame, typename ThisFrame>
SpecificEnergy BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::
GravitationalPotential(Instant const& t,
                       Position<InertialFrame> const& q) const {
  return ephemeris_->ComputeGravitationalPotential(q, t);
}

template<typename InertialFrame, typename ThisFrame>
AcceleratedRigidMotion<InertialFrame, ThisFrame>
BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::MotionOfThisFrame(
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

}  // namespace internal
}  // namespace _body_surface_reference_frame
}  // namespace physics
}  // namespace principia
