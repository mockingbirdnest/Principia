
// The files containing the tree of child classes of |DynamicFrame| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
#include "physics/dynamic_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_BARYCENTRIC_ROTATING_DYNAMIC_FRAME_HPP_
#define PRINCIPIA_PHYSICS_BARYCENTRIC_ROTATING_DYNAMIC_FRAME_HPP_

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
namespace internal_barycentric_rotating_dynamic_frame {

using base::not_null;
using geometry::AngularVelocity;
using geometry::Instant;
using geometry::Rotation;
using geometry::Vector;
using quantities::Acceleration;

template<typename InertialFrame, typename ThisFrame>
class BarycentricRotatingDynamicFrame
    : public DynamicFrame<InertialFrame, ThisFrame> {
 public:
  BarycentricRotatingDynamicFrame(
      not_null<Ephemeris<InertialFrame> const*> const ephemeris,
      not_null<MassiveBody const*> const primary,
      not_null<MassiveBody const*> const secondary);

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;
  RigidMotion<ThisFrame, InertialFrame> FromThisFrameAtTime(
      Instant const& t) const override;
  Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const override;

  void WriteToMessage(
      not_null<serialization::DynamicFrame*> const message) const override;

  static not_null<std::unique_ptr<BarycentricRotatingDynamicFrame>>
      ReadFromMessage(
          not_null<Ephemeris<InertialFrame> const*> const ephemeris,
          serialization::BarycentricRotatingDynamicFrame const& message);

 private:
  // Fills |*rotation| with the rotation that maps the basis of |InertialFrame|
  // to the basis of |ThisFrame|.  Fills |*angular_frequency| with the
  // corresponding angular velocity.
  static void ComputeAngularDegreesOfFreedom(
      DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
      DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
      not_null<Rotation<InertialFrame, ThisFrame>*> const rotation,
      not_null<AngularVelocity<InertialFrame>*> const angular_velocity);

  not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
  not_null<MassiveBody const*> const primary_;
  not_null<MassiveBody const*> const secondary_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const
      primary_trajectory_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const
      secondary_trajectory_;
  mutable typename ContinuousTrajectory<InertialFrame>::Hint primary_hint_;
  mutable typename ContinuousTrajectory<InertialFrame>::Hint secondary_hint_;
};

}  // namespace internal_barycentric_rotating_dynamic_frame

using internal_barycentric_rotating_dynamic_frame::
    BarycentricRotatingDynamicFrame;

}  // namespace physics
}  // namespace principia

#include "physics/barycentric_rotating_dynamic_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BARYCENTRIC_ROTATING_DYNAMIC_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
