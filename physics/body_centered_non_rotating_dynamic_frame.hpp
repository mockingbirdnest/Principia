#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using geometry::Instant;
using geometry::Vector;
using quantities::Acceleration;

namespace physics {

template<typename InertialFrame, typename ThisFrame>
class BodyCentredNonRotatingDynamicFrame
    : public DynamicFrame<InertialFrame, ThisFrame> {
 public:
  explicit BodyCentredNonRotatingDynamicFrame(
      ContinuousTrajectory<InertialFrame> const& centre_trajectory);

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;
  RigidMotion<ThisFrame, InertialFrame> FromThisFrameAtTime(
      Instant const& t) const override;
  Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const override;
};

}  // namespace physics
}  // namespace principia
