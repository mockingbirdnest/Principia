
#ifndef PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
#define PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_

#include "geometry/frame.hpp"
#include "geometry/rotation.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_dynamic_frame {

using base::not_null;
using geometry::Instant;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using quantities::Acceleration;

// The Frenet frame of a free fall trajectory in |Frame|.
// TODO(egg): this should actually depend on its template parameter somehow.
template<typename Frame>
using Frenet = geometry::Frame<serialization::Frame::PhysicsTag,
                               serialization::Frame::FRENET,
                               /*frame_is_inertial=*/false>;

// The definition of a reference frame |ThisFrame| in arbitrary motion with
// respect to the inertial reference frame |InertialFrame|.
template<typename InertialFrame, typename ThisFrame>
class DynamicFrame {
  static_assert(InertialFrame::is_inertial, "InertialFrame must be inertial");

 public:
  virtual ~DynamicFrame() = default;

  // At least one of |ToThisFrameAtTime| and |FromThisFrameAtTime| must be
  // overriden in derived classes; the default implementation inverts the other
  // one.
  virtual RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const;
  virtual RigidMotion<ThisFrame, InertialFrame> FromThisFrameAtTime(
      Instant const& t) const;

  // The acceleration due to the non-inertial motion of |ThisFrame| and gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is |GeometricAcceleration|.
  virtual Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const;

  // The definition of the Frenet frame of a free fall trajectory in |ThisFrame|
  // with the given |degrees_of_freedom| at instant |t|.
  virtual Rotation<Frenet<ThisFrame>, ThisFrame> FrenetFrame(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const;

  virtual void WriteToMessage(
      not_null<serialization::DynamicFrame*> const message) const = 0;

  // Dispatches to one of the subclasses depending on the contents of the
  // message.  Returns |nullptr| if no dynamic frame extension is found.
  static std::unique_ptr<DynamicFrame>
      ReadFromMessage(not_null<Ephemeris<InertialFrame> const*> ephemeris,
                      serialization::DynamicFrame const& message);

 private:
  virtual Vector<Acceleration, InertialFrame> GravitationalAcceleration(
      Instant const& t,
      Position<InertialFrame> const& q) const = 0;
  virtual AcceleratedRigidMotion<InertialFrame, ThisFrame> MotionOfThisFrame(
      Instant const& t) const = 0;
};

}  // namespace internal_dynamic_frame

using internal_dynamic_frame::DynamicFrame;
using internal_dynamic_frame::Frenet;

}  // namespace physics
}  // namespace principia

#include "physics/dynamic_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
