// The files containing the tree of child classes of |ReferenceFrame| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
#include "physics/reference_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
#define PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "physics/ephemeris.hpp"
#include "physics/reference_frame.hpp"
#include "physics/rigid_motion.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace _rigid_reference_frame {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_similar_motion;
using namespace principia::quantities::_named_quantities;

// The definition of a reference frame |ThisFrame| in arbitrary motion with
// respect to the inertial reference frame |InertialFrame|.
template<typename InertialFrame, typename ThisFrame>
class RigidReferenceFrame : public ReferenceFrame<InertialFrame, ThisFrame> {
 public:
  virtual ~RigidReferenceFrame() = default;

  SimilarMotion<InertialFrame, ThisFrame> ToThisFrameAtTimeSimilarly(
      Instant const& t) const final;
  SimilarMotion<ThisFrame, InertialFrame> FromThisFrameAtTimeSimilarly(
      Instant const& t) const final;

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

  // Dispatches to one of the subclasses depending on the contents of the
  // message.
  static not_null<std::unique_ptr<RigidReferenceFrame>>
      ReadFromMessage(serialization::ReferenceFrame const& message,
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

  // A helper function for computing the rotational movement of a frame defined
  // by two bodies.  Fills |rotation| with the rotation that maps the basis of
  // |InertialFrame| to the basis of |ThisFrame|.  Fills |angular_velocity| with
  // the corresponding angular velocity.
  static void ComputeAngularDegreesOfFreedom(
      DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
      DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
      Vector<Acceleration, InertialFrame> const& primary_acceleration,
      Vector<Acceleration, InertialFrame> const& secondary_acceleration,
      Rotation<InertialFrame, ThisFrame>& rotation,
      AngularVelocity<InertialFrame>& angular_velocity);
};

}  // namespace internal

using internal::RigidReferenceFrame;

}  // namespace _rigid_reference_frame
}  // namespace physics
}  // namespace principia

#include "physics/rigid_reference_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
