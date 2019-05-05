
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
namespace physics {
namespace internal_body_centred_non_rotating_dynamic_frame {

using base::not_null;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Vector;
using quantities::Acceleration;

// The origin of the frame is the centre of mass of the body.  The Y axis is at
// the intersection of the equator and the XY plane of |InertialFrame|, in the
// direction 90° + α from the X axis of |InertialFrame|, where α is the right
// ascension of the |polar_axis|.  The Z axis is the |polar_axis|.  The X axis
// in on the equator so that the frame has the same orientation as
// |InertialFrame|.
// For a non-rotating body, the axes are the same as those of |InertialFrame|.
// With respect to figure 1 of the 2015 report of the IAU WGCCRE if |polar_axis|
// is the north pole, or figure 2 if |polar_axis| is the positive pole,
// - the Y axis is the node Q from the figure;
// - the Z axis is the pole Z from the figure.
// Note that, when |InertialFrame| is the ICRS, these axes make |ThisFrame| the
// usual celestial reference frame (the GCRS) if |centre| is the Earth, as the X
// axis will point towards ♈; however, if |centre| is another planet,
// |ThisFrame| will not have the axes of its natural celestial reference frame,
// as the X axis will not point towards its equinox.
// REMOVE BEFORE FLIGHT: Do we really want that? This seems like a very
// contrived way get something that only makes sense for the Earth with respect
// to the ICRS; in particular, if the XY plane of |InertialFrame| is the
// invariable plane, Q should naturally be the X axis, not the Y axis.
template<typename InertialFrame, typename ThisFrame>
class BodyCentredNonRotatingDynamicFrame
    : public DynamicFrame<InertialFrame, ThisFrame> {
 public:
  BodyCentredNonRotatingDynamicFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      not_null<MassiveBody const*> centre);

  not_null<MassiveBody const*> centre() const;

  Instant t_min() const override;
  Instant t_max() const override;

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::DynamicFrame*> message) const override;

  static not_null<std::unique_ptr<BodyCentredNonRotatingDynamicFrame>>
  ReadFromMessage(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      serialization::BodyCentredNonRotatingDynamicFrame const& message);

 private:
  Vector<Acceleration, InertialFrame> GravitationalAcceleration(
      Instant const& t,
      Position<InertialFrame> const& q) const override;
  AcceleratedRigidMotion<InertialFrame, ThisFrame> MotionOfThisFrame(
      Instant const& t) const override;

  not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
  not_null<MassiveBody const*> const centre_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const centre_trajectory_;
  OrthogonalMap<InertialFrame, ThisFrame> const orthogonal_map_;
};


}  // namespace internal_body_centred_non_rotating_dynamic_frame

using internal_body_centred_non_rotating_dynamic_frame::
    BodyCentredNonRotatingDynamicFrame;

}  // namespace physics
}  // namespace principia

#include "physics/body_centred_non_rotating_dynamic_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BODY_CENTRED_NON_ROTATING_DYNAMIC_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
