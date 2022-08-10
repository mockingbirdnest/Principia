
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
using geometry::Arbitrary;
using geometry::Handedness;
using geometry::Instant;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using quantities::Acceleration;
using quantities::SpecificEnergy;

// The Frenet frame of a free fall trajectory in |Frame|.
// TODO(egg): this should actually depend on its template parameter somehow.
template<typename Frame>
using Frenet = geometry::Frame<serialization::Frame::PhysicsTag,
                               Arbitrary,
                               Handedness::Right,
                               serialization::Frame::FRENET>;

// The definition of a reference frame |ThisFrame| in arbitrary motion with
// respect to the inertial reference frame |InertialFrame|.
template<typename InertialFrame, typename ThisFrame>
class DynamicFrame {
  static_assert(InertialFrame::is_inertial, "InertialFrame must be inertial");

 public:
  virtual ~DynamicFrame() = default;

  // The operations that take an |Instant| are valid in the range
  // [t_min, t_max].
  virtual Instant t_min() const = 0;
  virtual Instant t_max() const = 0;

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

  // The acceleration of a particle at rest in |ThisFrame| at the given
  // |position| owing to non-inertial motion of |ThisFrame| and gravity,
  // excluding components with a rotation.
  // Let r be the radial vector (from the origin of |ThisFrame|) corresponding
  // to |position|.
  // Let r ↦ r″ be the vector field of free fall accelerations from rest in
  // |ThisFrame| at t. This function returns
  //   a = r″ - (rot r″) r / 2.
  // In a rotating reference frame, this may equivalently be expressed using
  // the second derivative of position with respect to the parametrization on
  // the angle θ rather than time, which eliminates the Euler acceleration:
  //   a = θ′² d²r/dθ², starting from a rest defined as dr/dθ = 0.
  // Either way, the vector field a derives from a potential.
  virtual Vector<Acceleration, ThisFrame>
  RotationFreeGeometricAccelerationAtRest(
      Instant const& t,
      Position<ThisFrame> const& position) const;

  // Computes the (scalar) potential from which the acceleration given by
  // |RotationFreeGeometricAccelerationAtRest| derives.
  virtual SpecificEnergy GeometricPotential(
      Instant const& t,
      Position<ThisFrame> const& position) const;

  // The definition of the Frenet frame of a free fall trajectory in |ThisFrame|
  // with the given |degrees_of_freedom| at instant |t|.
  virtual Rotation<Frenet<ThisFrame>, ThisFrame> FrenetFrame(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const;

  virtual void WriteToMessage(
      not_null<serialization::DynamicFrame*> message) const = 0;

  // Dispatches to one of the subclasses depending on the contents of the
  // message.
  static not_null<std::unique_ptr<DynamicFrame>>
      ReadFromMessage(serialization::DynamicFrame const& message,
                      not_null<Ephemeris<InertialFrame> const*> ephemeris);

 private:
  void ComputeGeometricAccelerations(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom,
      Vector<Acceleration, ThisFrame>& gravitational_acceleration,
      Vector<Acceleration, ThisFrame>& linear_acceleration,
      Vector<Acceleration, ThisFrame>& coriolis_acceleration,
      Vector<Acceleration, ThisFrame>& centrifugal_acceleration,
      Vector<Acceleration, ThisFrame>& euler_acceleration) const;

  virtual Vector<Acceleration, InertialFrame> GravitationalAcceleration(
      Instant const& t,
      Position<InertialFrame> const& q) const = 0;
  virtual SpecificEnergy GravitationalPotential(
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
