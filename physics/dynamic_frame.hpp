#pragma once

#include "geometry/frame.hpp"
#include "geometry/rotation.hpp"
#include "physics/rigid_motion.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using geometry::Rotation;

namespace physics {

// The Frenet frame of a free fall trajectory in |Frame|.
// TODO(egg): this should actually depend on its template parameter somehow.
template<typename Frame>
using Frenet = geometry::Frame<serialization::Frame::PhysicsTag,
                                    serialization::Frame::FRENET,
                                    false /*frame_is_inertial*/>;

// The definition of a reference frame |ThisFrame| in arbitrary motion with
// respect to the inertial reference frame |InertialFrame|.
template<typename InertialFrame, typename ThisFrame>
class DynamicFrame {
  static_assert(InertialFrame::is_inertial, "InertialFrame must be inertial");

 public:
  virtual ~DynamicFrame() = default;
  virtual RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const = 0;
  virtual RigidMotion<ThisFrame, InertialFrame> FromThisFrameAtTime(
      Instant const& t) const = 0;

  // The acceleration due to the non-inertial motion of |Frame| and gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is |GeometricAcceleration|.
  virtual Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const = 0;

  // The definition of the Frenet frame of a free fall trajectory in |ThisFrame|
  // with the given |DegreesOfFreedom| at the given |Instant|.  Unless
  // specialized otherwise by child classes, the member functions of the result
  // fail when called at instants other than |t|.
  virtual Rotation<Frenet<ThisFrame>, ThisFrame> FrenetFrame(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const;
};

// An inertial frame.
template<typename OtherFrame, typename ThisFrame>
class InertialFrame : public DynamicFrame<OtherFrame, ThisFrame> {
 public:
  InertialFrame(Velocity<OtherFrame> const& velocity,
                Position<OtherFrame> const& origin_at_epoch,
                Instant const& epoch,
                OrthogonalMap<OtherFrame, ThisFrame> const& orthogonal_map,
                std::function<Vector<Acceleration, OtherFrame>(
                    Instant const& t,
                    Position<OtherFrame> const& q)> gravity);

  RigidMotion<OtherFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;
  RigidMotion<ThisFrame, OtherFrame> FromThisFrameAtTime(
      Instant const& t) const override;

  // The acceleration due to the non-inertial motion of |Frame| and gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is |GeometricAcceleration|.
  Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const override;
 private:
  Velocity<OtherFrame> const velocity_;
  Position<OtherFrame> const origin_at_epoch_;
  Instant const epoch_;
  OrthogonalMap<OtherFrame, ThisFrame> const orthogonal_map_;
  std::function<Vector<Acceleration, OtherFrame>(
      Instant const& t,
      Position<OtherFrame> const& q)> gravity_;
};

}  // namespace physics
}  // namespace principia

#include "physics/dynamic_frame_body.hpp"
