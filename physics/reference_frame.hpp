#ifndef PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
#define PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_

#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "physics/ephemeris.hpp"
#include "physics/similar_motion.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace _reference_frame {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::physics::_similar_motion;
using namespace principia::quantities::_named_quantities;

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
class ReferenceFrame {
  static_assert(InertialFrame::is_inertial, "InertialFrame must be inertial");

 public:
  virtual ~ReferenceFrame() = default;

  // The operations that take an |Instant| are valid in the range
  // [t_min, t_max].
  virtual Instant t_min() const = 0;
  virtual Instant t_max() const = 0;

  // At least one of |ToThisFrameAtTimeSimilarly| and
  // |FromThisFrameAtTimeSimilarly| must be overriden in derived classes; the
  // default implementation inverts the other one.
  virtual SimilarMotion<InertialFrame, ThisFrame> ToThisFrameAtTimeSimilarly(
      Instant const& t) const;
  virtual SimilarMotion<ThisFrame, InertialFrame> FromThisFrameAtTimeSimilarly(
      Instant const& t) const;

  // The acceleration due to the non-inertial motion of |ThisFrame| and gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is |GeometricAcceleration|.
  virtual Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const = 0;

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
      Position<ThisFrame> const& position) const = 0;

  // Computes the (scalar) potential from which the acceleration given by
  // |RotationFreeGeometricAccelerationAtRest| derives.
  virtual SpecificEnergy GeometricPotential(
      Instant const& t,
      Position<ThisFrame> const& position) const = 0;

  // The definition of the Frenet frame of a free fall trajectory in |ThisFrame|
  // with the given |degrees_of_freedom| at instant |t|.
  virtual Rotation<Frenet<ThisFrame>, ThisFrame> FrenetFrame(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const;

  virtual void WriteToMessage(
      not_null<serialization::ReferenceFrame*> message) const = 0;

  // Dispatches to one of the subclasses depending on the contents of the
  // message.
  static not_null<std::unique_ptr<ReferenceFrame>>
      ReadFromMessage(serialization::ReferenceFrame const& message,
                      not_null<Ephemeris<InertialFrame> const*> ephemeris);
};

}  // namespace internal

using internal::ReferenceFrame;
using internal::Frenet;

}  // namespace _reference_frame
}  // namespace physics
}  // namespace principia

#include "physics/reference_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
