
#pragma once

#include "physics/body_centred_non_rotating_dynamic_frame.hpp"

#include "geometry/identity.hpp"
#include "geometry/rotation.hpp"

namespace principia {
namespace physics {
namespace internal_body_centred_non_rotating_dynamic_frame {

using geometry::AngularVelocity;
using geometry::Identity;
using geometry::OrthogonalMap;
using geometry::Rotation;

template<typename InertialFrame, typename ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
BodyCentredNonRotatingDynamicFrame(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    not_null<MassiveBody const*> const centre)
    : ephemeris_(ephemeris),
      centre_(centre),
      centre_trajectory_(ephemeris_->trajectory(centre_)),
      orthogonal_map_(
          [this]() {
            // Note that we cannot do this by making |equatorial| and
            // |biequatorial| virtual members of |MassiveBody|, because that
            // class is not templatized on |InertialFrame|.
            auto const rotating_body =
                dynamic_cast<RotatingBody<InertialFrame> const*>(&*centre_);
            if (rotating_body == nullptr) {
              return Identity<InertialFrame, ThisFrame>().Forget();
            }
            // In coordinates, the third parameter is |polar_axis|, but we seem
            // to be a bit confused as to which of these things should be
            // vectors or bivectors here.
            // TODO(egg): Figure that out.
            return Rotation<InertialFrame, ThisFrame>(
                rotating_body->equatorial(),
                rotating_body->biequatorial(),
                Wedge(rotating_body->equatorial(),
                      rotating_body->biequatorial())).Forget();
          }()) {}

template<typename InertialFrame, typename ThisFrame>
not_null<MassiveBody const*>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::centre() const {
  return centre_;
}

template<typename InertialFrame, typename ThisFrame>
Instant BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::t_min()
    const {
  return centre_trajectory_->t_min();
}

template<typename InertialFrame, typename ThisFrame>
Instant BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::t_max()
    const {
  return centre_trajectory_->t_max();
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const centre_degrees_of_freedom =
      centre_trajectory_->EvaluateDegreesOfFreedom(t);

  RigidTransformation<InertialFrame, ThisFrame> const
      rigid_transformation(centre_degrees_of_freedom.position(),
                           ThisFrame::origin,
                           orthogonal_map_);
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             AngularVelocity<InertialFrame>(),
             centre_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
void BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
WriteToMessage(not_null<serialization::DynamicFrame*> const message) const {
  message->MutableExtension(
      serialization::BodyCentredNonRotatingDynamicFrame::extension)->set_centre(
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

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, InertialFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
    GravitationalAcceleration(Instant const& t,
                              Position<InertialFrame> const& q) const {
  return ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(q, t);
}

template<typename InertialFrame, typename ThisFrame>
AcceleratedRigidMotion<InertialFrame, ThisFrame>
BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::MotionOfThisFrame(
    Instant const& t) const {
  return AcceleratedRigidMotion<InertialFrame, ThisFrame>(
             ToThisFrameAtTime(t),
             /*angular_acceleration_of_to_frame=*/{},
             /*acceleration_of_to_frame_origin=*/ephemeris_->
                 ComputeGravitationalAccelerationOnMassiveBody(centre_, t));
}

}  // namespace internal_body_centred_non_rotating_dynamic_frame
}  // namespace physics
}  // namespace principia
