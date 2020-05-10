
#pragma once

#include <functional>
#include <type_traits>

#include "base/not_null.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rotation.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_rigid_motion {

using base::not_null;
using geometry::AffineMap;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::RigidTransformation;
using geometry::Vector;
using geometry::Velocity;
using quantities::Acceleration;
using quantities::Length;
using quantities::Variation;
using quantities::si::Radian;

// The instantaneous motion of |ToFrame| with respect to |FromFrame|.
// This is the derivative of a |RigidTransformation<FromFrame, ToFrame>|.
// In order to invert, the |RigidTransformation| is needed, and we need its
// linear part anyway, so we store it (and we forward its action on positions).
template<typename FromFrame, typename ToFrame>
class RigidMotion final {
 public:
  RigidMotion(
      RigidTransformation<FromFrame, ToFrame> rigid_transformation,
      AngularVelocity<FromFrame> const& angular_velocity_of_to_frame,
      Velocity<FromFrame> const& velocity_of_to_frame_origin);

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same_v<F, T>>>
  RigidMotion(
      RigidTransformation<FromFrame, ToFrame> const& rigid_transformation,
      AngularVelocity<ToFrame> const& angular_velocity_of_from_frame,
      Velocity<ToFrame> const& velocity_of_from_frame_origin);

  RigidTransformation<FromFrame, ToFrame> const& rigid_transformation() const;
  // Returns |rigid_transformation().linear_map()|.
  OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map() const;
  AngularVelocity<FromFrame> const& angular_velocity_of_to_frame() const;
  Velocity<FromFrame> const& velocity_of_to_frame_origin() const;

  DegreesOfFreedom<ToFrame> operator()(
      DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const;

  RigidMotion<ToFrame, FromFrame> Inverse() const;

  // A rigid motion expressing that |FromFrame| and |ToFrame| have the same
  // axes, origin, and instantaneous motion.
  // This function is enabled only if both frames have the same handedness (this
  // is a requirement of OrthogonalMap::Identity) and if the |motion| of
  // FromFrame is a special case of that of |ToFrame| (see the comments on
  // |FrameMotion|).
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<(F::handedness == T::handedness &&
                                        F::motion <= T::motion)>>
  static RigidMotion Identity();

  // A factory that construct a non-rotating motion using the given degrees of
  // freedom.  Useful e.g. for save compatibility.
  static RigidMotion MakeNonRotatingMotion(
      DegreesOfFreedom<ToFrame> const& degrees_of_freedom_of_from_frame_origin);

  void WriteToMessage(not_null<serialization::RigidMotion*> message) const;
  static RigidMotion ReadFromMessage(serialization::RigidMotion const& message);

 private:
  RigidTransformation<FromFrame, ToFrame> rigid_transformation_;
  // d/dt rigid_transformation⁻¹(basis of ToFrame). The positively oriented
  // orthogonal bases of |FromFrame| are acted upon faithfully and transitively
  // by SO(FromFrame), so this lies in the tangent space, i.e., the Lie algebra
  // 𝖘𝔬(FromFrame) ≅ FromFrame ∧ FromFrame.
  AngularVelocity<FromFrame> angular_velocity_of_to_frame_;
  // d/dt rigid_transformation⁻¹(ToFrame::origin).
  Velocity<FromFrame> velocity_of_to_frame_origin_;

  template<typename From, typename Through, typename To>
  friend RigidMotion<From, To> operator*(
      RigidMotion<Through, To> const& left,
      RigidMotion<From, Through> const& right);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame> operator*(
    RigidMotion<ThroughFrame, ToFrame> const& left,
    RigidMotion<FromFrame, ThroughFrame> const& right);

// A |RigidTransformation|, its first derivative (a |RigidMotion|), and its
// second derivative (angular and linear accelerations).
template<typename FromFrame, typename ToFrame>
class AcceleratedRigidMotion final {
 public:
  AcceleratedRigidMotion(
      RigidMotion<FromFrame, ToFrame> rigid_motion,
      Variation<AngularVelocity<FromFrame>> const&
          angular_acceleration_of_to_frame,
      Vector<Acceleration, FromFrame> const& acceleration_of_to_frame_origin);

  RigidMotion<FromFrame, ToFrame> const& rigid_motion() const;
  Variation<AngularVelocity<FromFrame>> const&
  angular_acceleration_of_to_frame() const;
  Vector<Acceleration, FromFrame> const& acceleration_of_to_frame_origin()
      const;

 private:
  RigidMotion<FromFrame, ToFrame> const rigid_motion_;
  // d/dt rigid_motion_.angular_velocity_of_to_frame().
  Variation<AngularVelocity<FromFrame>> const angular_acceleration_of_to_frame_;
  // d/dt rigid_motion_.velocity_of_to_frame_origin().
  Vector<Acceleration, FromFrame> const acceleration_of_to_frame_origin_;
};

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         RigidMotion<FromFrame, ToFrame> const& rigid_motion);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    AcceleratedRigidMotion<FromFrame, ToFrame> const& accelerated_rigid_motion);

}  // namespace internal_rigid_motion

using internal_rigid_motion::AcceleratedRigidMotion;
using internal_rigid_motion::RigidMotion;
using internal_rigid_motion::RigidTransformation;

}  // namespace physics
}  // namespace principia

#include "physics/rigid_motion_body.hpp"
