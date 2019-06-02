
// The files containing the tree of child classes of |DynamicFrame| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
#include "physics/dynamic_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_BODY_SURFACE_DYNAMIC_FRAME_HPP_
#define PRINCIPIA_PHYSICS_BODY_SURFACE_DYNAMIC_FRAME_HPP_

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rotating_body.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_body_surface_dynamic_frame {

using base::not_null;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using quantities::Acceleration;

// The origin of the frame is the centre of mass of the body.  The X axis is at
// the intersection of the equator and the prime meridian.  The Z axis is the
// |polar_axis|.  The Y axis in on the equator so that the frame has the same
// orientation as |InertialFrame|.
// The X, Y, Z axes are the same as those shown on figure 1 of the 2015 report
// of the IAU WGCCRE if |polar_axis| is the north pole, or figure 2 if
// |polar_axis| is the positive pole.
template<typename InertialFrame, typename ThisFrame>
class BodySurfaceDynamicFrame
    : public DynamicFrame<InertialFrame, ThisFrame> {
 public:
  BodySurfaceDynamicFrame(not_null<Ephemeris<InertialFrame> const*> ephemeris,
                          not_null<RotatingBody<InertialFrame> const*> centre);

  not_null<RotatingBody<InertialFrame> const*> centre() const;

  Instant t_min() const override;
  Instant t_max() const override;

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::DynamicFrame*> message) const override;

  static not_null<std::unique_ptr<BodySurfaceDynamicFrame>> ReadFromMessage(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      serialization::BodySurfaceDynamicFrame const& message);

 private:
  Vector<Acceleration, InertialFrame> GravitationalAcceleration(
      Instant const& t,
      Position<InertialFrame> const& q) const override;
  AcceleratedRigidMotion<InertialFrame, ThisFrame> MotionOfThisFrame(
      Instant const& t) const override;

  not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
  not_null<RotatingBody<InertialFrame> const*> const centre_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const centre_trajectory_;
};

}  // namespace internal_body_surface_dynamic_frame

using internal_body_surface_dynamic_frame::BodySurfaceDynamicFrame;

}  // namespace physics
}  // namespace principia

#include "physics/body_surface_dynamic_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BODY_SURFACE_DYNAMIC_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
