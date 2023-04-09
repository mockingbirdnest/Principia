// The files containing the tree of child classes of |RigidReferenceFrame| must
// be included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
#include "physics/rigid_reference_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_BODY_SURFACE_REFERENCE_FRAME_HPP_
#define PRINCIPIA_PHYSICS_BODY_SURFACE_REFERENCE_FRAME_HPP_

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rotating_body.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace _body_surface_reference_frame {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::quantities::_named_quantities;

// The origin of the frame is the centre of mass of the body.  The X axis is at
// the intersection of the equator and the prime meridian.  The Z axis is the
// |polar_axis|.  The Y axis in on the equator so that the frame has the same
// orientation as |InertialFrame|.
// The X, Y, Z axes are the same as those shown on figure 1 of the 2015 report
// of the IAU WGCCRE if |polar_axis| is the north pole, or figure 2 if
// |polar_axis| is the positive pole.
template<typename InertialFrame, typename ThisFrame>
class BodySurfaceReferenceFrame : public RigidReferenceFrame<InertialFrame,
                                                             ThisFrame> {
  static_assert(ThisFrame::may_rotate);

 public:
  BodySurfaceReferenceFrame(not_null<Ephemeris<InertialFrame> const*> ephemeris,
                          not_null<RotatingBody<InertialFrame> const*> centre);

  not_null<RotatingBody<InertialFrame> const*> centre() const;

  Instant t_min() const override;
  Instant t_max() const override;

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::ReferenceFrame*> message) const override;

  static not_null<std::unique_ptr<BodySurfaceReferenceFrame>> ReadFromMessage(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      serialization::BodySurfaceReferenceFrame const& message);

 private:
  Vector<Acceleration, InertialFrame> GravitationalAcceleration(
      Instant const& t,
      Position<InertialFrame> const& q) const override;
  SpecificEnergy GravitationalPotential(
      Instant const& t,
      Position<InertialFrame> const& q) const override;
  AcceleratedRigidMotion<InertialFrame, ThisFrame> MotionOfThisFrame(
      Instant const& t) const override;

  not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
  not_null<RotatingBody<InertialFrame> const*> const centre_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const centre_trajectory_;
};

}  // namespace internal

using internal::BodySurfaceReferenceFrame;

}  // namespace _body_surface_reference_frame
}  // namespace physics
}  // namespace principia

namespace principia::physics {
using namespace principia::physics::_body_surface_reference_frame;
}  // namespace principia::physics

#include "physics/body_surface_reference_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BODY_SURFACE_REFERENCE_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
