#pragma once

#include <functional>

#include "geometry/affine_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rotation.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using geometry::AffineMap;
using geometry::AngularVelocity;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Vector;
using quantities::Acceleration;
using quantities::Length;
using quantities::si::Radian;

namespace physics {

// An arbitrary rigid transformation.  Simultaneous positions between two frames
// are always related by such a transformation.
template<typename FromFrame, typename ToFrame>
using RigidTransformation =
    AffineMap<FromFrame, ToFrame, Length, OrthogonalMap>;

// The instantaneous motion of |ToFrame| with respect to |FromFrame|.
// This is the derivative of a |RigidTransformation<FromFrame, ToFrame>|.
// In order to invert, the |RigidTransformation| is needed, and we need its
// linear part anyway, so we store it (and we forward its action on positions).
template<typename FromFrame, typename ToFrame>
class RigidMotion {
 public:
  RigidMotion(
      RigidTransformation<FromFrame, ToFrame> const& rigid_transformation,
      AngularVelocity<FromFrame> const& rotation,
      Velocity<FromFrame> const& translation);
  RigidMotion(
      RigidTransformation<FromFrame, ToFrame> const& rigid_transformation,
      AngularVelocity<ToFrame> const& rotation,
      Velocity<ToFrame> const& translation);
  ~RigidMotion() = default;

  RigidTransformation<FromFrame, ToFrame> const& rigid_transformation() const;
  // Returns |rigid_transformation().linear_map()|.
  OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map() const;

  DegreesOfFreedom<ToFrame> operator()(
      DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const;

  RigidMotion<ToFrame, FromFrame> Inverse() const;

 private:
  RigidTransformation<FromFrame, ToFrame> const rigid_transformation_;
  // d/dt rigid_transformation(basis of FromFrame). The positively oriented
  // orthogonal bases of |ToFrame| are acted upon faithfully and transitively by
  // SO(ToFrame), so this lies in the tangent space, i.e., the Lie algebra
  // 𝖘𝔬(ToFrame) ≅ ToFrame ∧ ToFrame.
  AngularVelocity<ToFrame> rotation_;
  // d/dt rigid_transformation(FromFrame::origin).
  Velocity<ToFrame> translation_;

  template<typename From, typename Through, typename To>
  friend RigidMotion<From, To> operator*(
      RigidMotion<Through, To> const& left,
      RigidMotion<From, Through> const& right);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame> operator*(
    RigidMotion<ThroughFrame, ToFrame> const& left,
    RigidMotion<FromFrame, ThroughFrame> const& right);

}  // namespace physics
}  // namespace principia

#include "physics/rigid_motion_body.hpp"
