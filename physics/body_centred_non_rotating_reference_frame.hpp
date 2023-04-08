// The files containing the tree of child classes of |RigidReferenceFrame| must
// be included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
#include "physics/rigid_reference_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_BODY_CENTRED_NON_ROTATING_REFERENCE_FRAME_HPP_
#define PRINCIPIA_PHYSICS_BODY_CENTRED_NON_ROTATING_REFERENCE_FRAME_HPP_

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/space.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace _body_centred_non_rotating_reference_frame {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_space;
using namespace principia::quantities::_named_quantities;

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
// axis will point towards ♈︎; however, if |centre| is another planet,
// |ThisFrame| will not have the axes of its natural celestial reference frame,
// as the X axis will not point towards its equinox.
template<typename InertialFrame, typename ThisFrame>
class BodyCentredNonRotatingReferenceFrame
    : public RigidReferenceFrame<InertialFrame, ThisFrame> {
  static_assert(!ThisFrame::is_inertial);

 public:
  BodyCentredNonRotatingReferenceFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      not_null<MassiveBody const*> centre);

  not_null<MassiveBody const*> centre() const;

  Instant t_min() const override;
  Instant t_max() const override;

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::RigidReferenceFrame*> message) const override;

  static not_null<std::unique_ptr<BodyCentredNonRotatingReferenceFrame>>
  ReadFromMessage(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      serialization::BodyCentredNonRotatingReferenceFrame const& message);

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
  not_null<MassiveBody const*> const centre_;
  not_null<ContinuousTrajectory<InertialFrame> const*> const centre_trajectory_;
  OrthogonalMap<InertialFrame, ThisFrame> const orthogonal_map_;
};


}  // namespace internal

using internal::BodyCentredNonRotatingReferenceFrame;

}  // namespace _body_centred_non_rotating_reference_frame
}  // namespace physics
}  // namespace principia

namespace principia::physics {
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
}  // namespace principia::physics

#include "physics/body_centred_non_rotating_reference_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BODY_CENTRED_NON_ROTATING_REFERENCE_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
