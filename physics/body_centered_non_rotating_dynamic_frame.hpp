
// The files containing the tree of child classes of |DynamicFrame| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
#include "physics/dynamic_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_BODY_CENTRED_NON_ROTATING_DYNAMIC_FRAME_HPP_
#define PRINCIPIA_PHYSICS_BODY_CENTRED_NON_ROTATING_DYNAMIC_FRAME_HPP_

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
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
  BodyCentredNonRotatingDynamicFrame(
      not_null<Ephemeris<InertialFrame> const*> const ephemeris,
      not_null<MassiveBody const*> const centre);

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;
  RigidMotion<ThisFrame, InertialFrame> FromThisFrameAtTime(
      Instant const& t) const override;
  Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const override;

  void WriteToMessage(
      not_null<serialization::DynamicFrame*> const message) const override;

  static not_null<std::unique_ptr<BodyCentredNonRotatingDynamicFrame>>
      ReadFromMessage(
          not_null<Ephemeris<InertialFrame> const*> const ephemeris,
          serialization::BodyCentredNonRotatingDynamicFrame const& message);

 private:
  not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
  not_null<MassiveBody const*> const centre_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const centre_trajectory_;
  mutable typename ContinuousTrajectory<InertialFrame>::Hint hint_;
};

}  // namespace physics
}  // namespace principia

#include "physics/body_centered_non_rotating_dynamic_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BODY_CENTRED_NON_ROTATING_DYNAMIC_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
