#pragma once

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

using base::not_null;
using geometry::Bivector;
using geometry::Instant;
using geometry::Rotation;
using geometry::Vector;
using quantities::Acceleration;
using quantities::AngularFrequency;

namespace physics {

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

 private:
  // Fills |*rotation| with the rotation that maps the basis of |InertialFrame|
  // to the basis of |ThisFrame|.  Fills |*angular_frequency| with the
  // corresponding angular velocity.  |barycentre_degrees_of_freedom| must be a
  // convex combination of the two other degrees of freedom.
  static void FromBasisOfInertialFrameToBasisOfThisFrame(
      DegreesOfFreedom<InertialFrame> const& barycentre_degrees_of_freedom,
      DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
      DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
      not_null<Rotation<InertialFrame, ThisFrame>*> const rotation,
      not_null<Bivector<AngularFrequency, InertialFrame>*> const
          angular_frequency);

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

}  // namespace physics
}  // namespace principia

#include "physics/barycentric_rotating_dynamic_frame_body.hpp"
