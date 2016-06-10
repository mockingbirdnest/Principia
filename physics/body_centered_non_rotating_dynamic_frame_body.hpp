
#pragma once

#include "physics/body_centered_non_rotating_dynamic_frame.hpp"

#include "geometry/identity.hpp"

namespace principia {
namespace physics {
namespace internal_body_centred_non_rotating_dynamic_frame {

using geometry::AngularVelocity;
using geometry::Identity;

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
                           Identity<InertialFrame, ThisFrame>().Forget());
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             AngularVelocity<InertialFrame>(),
             centre_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<ThisFrame, InertialFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
FromThisFrameAtTime(Instant const& t) const {
  return ToThisFrameAtTime(t).Inverse();
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
GeometricAcceleration(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  auto const to_this_frame = ToThisFrameAtTime(t);
  auto const from_this_frame = to_this_frame.Inverse();

  Vector<Acceleration, ThisFrame> const gravitational_acceleration_at_point =
      to_this_frame.orthogonal_map()(
          ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(
              from_this_frame.rigid_transformation()(
                  degrees_of_freedom.position()), t));
  Vector<Acceleration, ThisFrame> const linear_acceleration =
      to_this_frame.orthogonal_map()(
          -ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(
              centre_, t));

  Vector<Acceleration, ThisFrame> const& fictitious_acceleration =
      linear_acceleration;
  return gravitational_acceleration_at_point + fictitious_acceleration;
}

template<typename InertialFrame, typename ThisFrame>
void BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
WriteToMessage(not_null<serialization::DynamicFrame*> const message) const {
  message->MutableExtension(
      serialization::BodyCentredNonRotatingDynamicFrame::
          body_centred_non_rotating_dynamic_frame)->set_centre(
              ephemeris_->serialization_index_for_body(centre_));
}

template<typename InertialFrame, typename ThisFrame>
not_null<std::unique_ptr<
    BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>>>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::BodyCentredNonRotatingDynamicFrame const& message) {
  return std::make_unique<BodyCentredNonRotatingDynamicFrame>(
             ephemeris,
             ephemeris->body_for_serialization_index(message.centre()));
}

}  // namespace internal_body_centred_non_rotating_dynamic_frame
}  // namespace physics
}  // namespace principia
