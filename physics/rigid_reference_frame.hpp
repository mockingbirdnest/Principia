// The files containing the tree of child classes of `ReferenceFrame` must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
#include "physics/reference_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
#define PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_

#include <memory>

#include "base/algebra.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/similar_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace _rigid_reference_frame {
namespace internal {

using namespace principia::base::_algebra;
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
using namespace principia::quantities::_quantities;

// The definition of a reference frame `ThisFrame` in arbitrary motion with
// respect to the inertial reference frame `InertialFrame`.
template<typename InertialFrame, typename ThisFrame>
class RigidReferenceFrame : public ReferenceFrame<InertialFrame, ThisFrame> {
 public:
  ~RigidReferenceFrame() override = default;

  SimilarMotion<InertialFrame, ThisFrame> ToThisFrameAtTimeSimilarly(
      Instant const& t) const final;
  SimilarMotion<ThisFrame, InertialFrame> FromThisFrameAtTimeSimilarly(
      Instant const& t) const final;

  // At least one of `ToThisFrameAtTime` and `FromThisFrameAtTime` must be
  // overriden in derived classes; the default implementation inverts the other
  // one.
  virtual RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const;
  virtual RigidMotion<ThisFrame, InertialFrame> FromThisFrameAtTime(
      Instant const& t) const;

  // The acceleration due to the non-inertial motion of `ThisFrame` and gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is `GeometricAcceleration`.
  Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const override;

  // The acceleration of a particle at rest in `ThisFrame` at the given
  // `position` owing to non-inertial motion of `ThisFrame` and gravity,
  // excluding components with a rotation.
  // Let r be the radial vector (from the origin of `ThisFrame`) corresponding
  // to `position`.
  // Let r ↦ r″ be the vector field of free fall accelerations from rest in
  // `ThisFrame` at t. This function returns
  //   a = r″ - (rot r″) r / 2.
  // In a rotating reference frame, this may equivalently be expressed using
  // the second derivative of position with respect to the parametrization on
  // the angle θ rather than time, which eliminates the Euler acceleration:
  //   a = θ′² d²r/dθ², starting from a rest defined as dr/dθ = 0.
  // Either way, the vector field a derives from a potential.
  Vector<Acceleration, ThisFrame>
  RotationFreeGeometricAccelerationAtRest(
      Instant const& t,
      Position<ThisFrame> const& position) const override;

  // Computes the (scalar) potential from which the acceleration given by
  // `RotationFreeGeometricAccelerationAtRest` derives.
  SpecificEnergy GeometricPotential(
      Instant const& t,
      Position<ThisFrame> const& position) const override;

  // Dispatches to one of the subclasses depending on the contents of the
  // message.
  static not_null<std::unique_ptr<RigidReferenceFrame>>
      ReadFromMessage(serialization::ReferenceFrame const& message,
                      not_null<Ephemeris<InertialFrame> const*> ephemeris);

 protected:
  // A helper function for computing the rotational movement of a frame defined
  // by two bodies.
  static void ComputeAngularDegreesOfFreedom(
      DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
      DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
      Vector<Acceleration, InertialFrame> const& primary_acceleration,
      Vector<Acceleration, InertialFrame> const& secondary_acceleration,
      Rotation<InertialFrame, ThisFrame>& rotation,
      AngularVelocity<InertialFrame>& angular_velocity);

  // The component names are somewhat notional.  This is not a Frenet frame and
  // anyway, the derivative trihedra are not even orthogonal.
  template<typename ScalarF, typename ScalarB, int order = 0>
  struct Trihedron {
    Vector<Derivative<ScalarF, Time, order>, InertialFrame> fore;
    Vector<Derivative<Product<ScalarF, ScalarB>, Time, order>,
           InertialFrame> normal;
    Bivector<Derivative<ScalarB, Time, order>, InertialFrame> binormal;
  };

  // Computes the orthogonal and orthonormal trihedra associated with
  // `ThisFrame`.
  static void ComputeTrihedra(
      Displacement<InertialFrame> const& r,
      Velocity<InertialFrame> const& ṙ,
      Trihedron<Length, ArealSpeed>& orthogonal,
      Trihedron<double, double>& orthonormal);

  // Computes the first derivative of the preceding trihedra.
  static void ComputeTrihedraDerivatives(
      Displacement<InertialFrame> const& r,
      Velocity<InertialFrame> const& ṙ,
      Vector<Acceleration, InertialFrame> const& r̈,
      Trihedron<Length, ArealSpeed> const& orthogonal,
      Trihedron<double, double> const& orthonormal,
      Trihedron<Length, ArealSpeed, 1>& 𝛛orthogonal,
      Trihedron<double, double, 1>& 𝛛orthonormal);

  // Computes the second derivative of the preceding trihedra.
  static void ComputeTrihedraDerivatives2(
      Displacement<InertialFrame> const& r,
      Velocity<InertialFrame> const& ṙ,
      Vector<Acceleration, InertialFrame> const& r̈,
      Vector<Jerk, InertialFrame> const& r⁽³⁾,
      Trihedron<Length, ArealSpeed> const& orthogonal,
      Trihedron<double, double> const& orthonormal,
      Trihedron<Length, ArealSpeed, 1> const& 𝛛orthogonal,
      Trihedron<double, double, 1> const& 𝛛orthonormal,
      Trihedron<Length, ArealSpeed, 2>& 𝛛²orthogonal,
      Trihedron<double, double, 2>& 𝛛²orthonormal);

  // Computes the rotation that maps `InertialFrame` to `ThisFrame`.
  static Rotation<InertialFrame, ThisFrame> ComputeRotation(
      Trihedron<double, double> const& orthonormal);

  // Computes the angular velocity of `ThisFrame` in `InertialFrame`.
  static AngularVelocity<InertialFrame> ComputeAngularVelocity(
      Trihedron<double, double> const& orthonormal,
      Trihedron<double, double, 1> const& 𝛛orthonormal);

  // Computes the angular acceleration of `ThisFrame` in `InertialFrame`.
  static Bivector<AngularAcceleration, InertialFrame>
  ComputeAngularAcceleration(Trihedron<double, double> const& orthonormal,
                             Trihedron<double, double, 1> const& 𝛛orthonormal,
                             Trihedron<double, double, 2> const& 𝛛²orthonormal);

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

}  // namespace internal

using internal::RigidReferenceFrame;

}  // namespace _rigid_reference_frame
}  // namespace physics
}  // namespace principia

#include "physics/rigid_reference_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
