#pragma once

#include <functional>

#include "geometry/affine_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using geometry::AffineMap;
using geometry::AngularVelocity;
using geometry::Instant;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using quantities::Acceleration;
using quantities::Length;

namespace physics {

// TODO(egg): AffineMap needs an accessor for its linear map.
// TODO(egg): why not OrthogonalMap?
template<typename FromFrame, typename ToFrame>
using RigidTransformation = AffineMap<FromFrame, ToFrame, Length, Rotation>;

// The instantaneous motion of |ToFrame| with respect to |FromFrame|.
// This is the derivative of a |RigidTransformation<FromFrame, ToFrame>|.
// In order to invert, the |RigidTransformation| is needed, and we need its
// linear part anyway, so we store it (and we forward its action on positions).
template<typename FromFrame, typename ToFrame>
class RigidMotion {
 public:
  RigidTransformation<FromFrame, ToFrame> const& rigid_transformation();
  // Returns |rigid_transformation().linear_map()|.
  Rotation<FromFrame, ToFrame> const& rotation() const;

  DegreesOfFreedom<ToFrame> operator()(
      DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const;

  RigidMotion<ToFrame, FromFrame> Inverse() const;

 private:
  RigidTransformation<FromFrame, ToFrame> rigid_tranformation_;
  // d/dt rigid_transformation?¹(basis of ToFrame). The positively oriented
  // orthogonal bases of FromFrame are acted upon faithfully and transitively by
  // |Rotation<FromFrame, FromFrame>|, so this lies in the tangent space, i.e.,
  // the Lie algebra.
  AngularVelocity<FromFrame> rotation_;
  // d/dt rigid_transformation?¹(ToFrame::origin).
  Velocity<FromFrame> translation_;
};

// The definition of a reference frame |Frame| in arbitrary motion with respect
// to |BaseFrame|.
template<typename BaseFrame, typename Frame>
class DynamicFrame {
 public:
  virtual RigidMotion<BaseFrame, Frame> ToFrameAtTime(
      Instant const& t) = 0;
  virtual RigidMotion<Frame, BaseFrame> FromFrameAtTime(
      Instant const& t) = 0;

  // The acceleration due to the non-inertial motion of |Frame| and gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is |GeometricAcceleration|.
  virtual Vector<Acceleration, Frame> GeometricAcceleration(
      Instant const& t,
      Position const& q) = 0;
};

}  // namespace physics
}  // namespace principia
