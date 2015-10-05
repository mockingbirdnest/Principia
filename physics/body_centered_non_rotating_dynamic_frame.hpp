#pragma once

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/massive_body.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using base::not_null;
using geometry::Instant;
using geometry::Vector;
using quantities::Acceleration;

namespace physics {

template<typename InertialFrame, typename ThisFrame>
class BodyCentredNonRotatingDynamicFrame
    : public DynamicFrame<InertialFrame, ThisFrame> {
 public:
  explicit BodyCentredNonRotatingDynamicFrame(
      not_null<Ephemeris<InertialFrame> const*> const ephemeris,
      not_null<MassiveBody const*> const centre);

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;
  RigidMotion<ThisFrame, InertialFrame> FromThisFrameAtTime(
      Instant const& t) const override;
  Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const override;

 private:
    not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
    not_null<MassiveBody const*> const centre_;
    not_null<ContinuousTrajectory<InertialFrame> const*> const
        centre_trajectory_;
    ContinuousTrajectory<InertialFrame> hint_;
};

}  // namespace physics
}  // namespace principia
