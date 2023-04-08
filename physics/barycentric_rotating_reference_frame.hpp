// The files containing the tree of child classes of |RigidReferenceFrame| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
#include "physics/rigid_reference_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_BARYCENTRIC_ROTATING_REFERENCE_FRAME_HPP_
#define PRINCIPIA_PHYSICS_BARYCENTRIC_ROTATING_REFERENCE_FRAME_HPP_

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace _barycentric_rotating_reference_frame {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::quantities::_named_quantities;

// The origin of the frame is the barycentre of the two bodies.  The X axis
// points to the secondary.  The Y axis is in the direction of the velocity of
// the secondary with respect to the primary.  The Z axis is in the direction of
// the angular velocity of the system.  The basis has the same orientation as
// |InertialFrame|.
template<typename InertialFrame, typename ThisFrame>
class BarycentricRotatingReferenceFrame
    : public RigidReferenceFrame<InertialFrame, ThisFrame> {
  static_assert(ThisFrame::may_rotate);

 public:
  BarycentricRotatingReferenceFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      not_null<MassiveBody const*> primary,
      not_null<MassiveBody const*> secondary);

  not_null<MassiveBody const*> primary() const;
  not_null<MassiveBody const*> secondary() const;

  Instant t_min() const override;
  Instant t_max() const override;

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::RigidReferenceFrame*> message) const override;

  static not_null<std::unique_ptr<BarycentricRotatingReferenceFrame>>
  ReadFromMessage(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      serialization::BarycentricRotatingReferenceFrame const& message);

 private:
  Vector<Acceleration, InertialFrame> GravitationalAcceleration(
      Instant const& t,
      Position<InertialFrame> const& q) const override;
  SpecificEnergy GravitationalPotential(
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
  not_null<MassiveBody const*> const primary_;
  not_null<MassiveBody const*> const secondary_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const
      primary_trajectory_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const
      secondary_trajectory_;
};

}  // namespace internal

using internal::BarycentricRotatingReferenceFrame;

}  // namespace _barycentric_rotating_reference_frame
}  // namespace physics
}  // namespace principia

namespace principia::physics {
using namespace principia::physics::_barycentric_rotating_reference_frame;
}  // namespace principia::physics

#include "physics/barycentric_rotating_reference_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BARYCENTRIC_ROTATING_REFERENCE_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
