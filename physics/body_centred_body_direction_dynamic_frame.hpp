
// The files containing the tree of child classes of |DynamicFrame| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
#include "physics/dynamic_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_BODY_CENTRED_BODY_DIRECTION_DYNAMIC_FRAME_HPP_
#define PRINCIPIA_PHYSICS_BODY_CENTRED_BODY_DIRECTION_DYNAMIC_FRAME_HPP_

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_body_centred_body_direction_dynamic_frame {

using base::not_null;
using geometry::AngularVelocity;
using geometry::Instant;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using quantities::Acceleration;

// The origin of the frame is the centre of mass of the primary body.  The X
// axis points to the secondary.  The Y axis is in the direction of the velocity
// of the secondary with respect to the primary.  The Z axis is in the direction
// of the angular velocity of the system.  The basis has the same orientation as
// |InertialFrame|.
template<typename InertialFrame, typename ThisFrame>
class BodyCentredBodyDirectionDynamicFrame
    : public DynamicFrame<InertialFrame, ThisFrame> {
  static_assert(ThisFrame::may_rotate);

 public:
  BodyCentredBodyDirectionDynamicFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      not_null<MassiveBody const*> primary,
      not_null<MassiveBody const*> secondary);

  BodyCentredBodyDirectionDynamicFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      std::function<Trajectory<InertialFrame> const&()> primary_trajectory,
      not_null<MassiveBody const*> secondary);

  not_null<MassiveBody const*> primary() const;
  not_null<MassiveBody const*> secondary() const;

  Instant t_min() const override;
  Instant t_max() const override;

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::DynamicFrame*> message) const override;

  static not_null<std::unique_ptr<BodyCentredBodyDirectionDynamicFrame>>
      ReadFromMessage(
          not_null<Ephemeris<InertialFrame> const*> ephemeris,
          serialization::BodyCentredBodyDirectionDynamicFrame const& message);

 private:
  Vector<Acceleration, InertialFrame> GravitationalAcceleration(
      Instant const& t,
      Position<InertialFrame> const& q) const override;
  AcceleratedRigidMotion<InertialFrame, ThisFrame> MotionOfThisFrame(
      Instant const& t) const override;

  // Fills |rotation| with the rotation that maps the basis of |InertialFrame|
  // to the basis of |ThisFrame|.  Fills |angular_velocity| with the
  // corresponding angular velocity.
  static void ComputeAngularDegreesOfFreedom(
      DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
      DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
      Rotation<InertialFrame, ThisFrame>& rotation,
      AngularVelocity<InertialFrame>& angular_velocity);

  not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
  MassiveBody const* const primary_;
  not_null<MassiveBody const*> const secondary_;
  std::function<Vector<Acceleration, InertialFrame>(
      Position<InertialFrame> const& position,
      Instant const& t)> compute_gravitational_acceleration_on_primary_;
  std::function<Trajectory<InertialFrame> const&()> const primary_trajectory_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const
      secondary_trajectory_;
};

}  // namespace internal_body_centred_body_direction_dynamic_frame

using internal_body_centred_body_direction_dynamic_frame::
    BodyCentredBodyDirectionDynamicFrame;

}  // namespace physics
}  // namespace principia

#include "physics/body_centred_body_direction_dynamic_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BODY_CENTRED_BODY_DIRECTION_DYNAMIC_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
