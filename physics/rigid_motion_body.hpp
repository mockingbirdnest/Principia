#pragma once

#include "physics/rigid_motion.hpp"

namespace principia {
namespace physics {

template<typename FromFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame>::RigidMotion(
    RigidTransformation<FromFrame, ToFrame> const& rigid_transformation,
    AngularVelocity<ToFrame> const& rotation,
    Velocity<ToFrame> const& translation)
    : rigid_transformation_(rigid_transformation),
      rotation_(rotation),
      translation_(translation) {}

template<typename FromFrame, typename ToFrame>
RigidTransformation<FromFrame, ToFrame> const&
RigidMotion<FromFrame, ToFrame>::rigid_transformation() const {
  return rigid_transformation_;
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> const&
RigidMotion<FromFrame, ToFrame>::orthogonal_map() const {
  return rigid_transformation_.linear_map();
}

template<typename FromFrame, typename ToFrame>
DegreesOfFreedom<ToFrame> RigidMotion<FromFrame, ToFrame>::operator()(
    DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const {
  return {rigid_transformation_(degrees_of_freedom.position()),
          orthogonal_map()(degrees_of_freedom.velocity()) + translation_ +
              rotation_ * orthogonal_map()(degrees_of_freedom.position() -
                                           FromFrame::origin) / Radian};
}

template<typename FromFrame, typename ToFrame>
RigidMotion<ToFrame, FromFrame>
RigidMotion<FromFrame, ToFrame>::Inverse() const {
  return RigidMotion<ToFrame, FromFrame>(
      rigid_transformation_.Inverse(),
      -orthogonal_map().Inverse()(rotation_),
      -orthogonal_map().Inverse()(translation_));
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame> operator*(
    RigidMotion<ThroughFrame, ToFrame> const& left,
    RigidMotion<FromFrame, ThroughFrame> const& right) {
  return RigidMotion<FromFrame, ToFrame>(
      left.rigid_transformation() * right.rigid_transformation(),
      left.rotation_ + left.orthogonal_map()(right.rotation_),
      left(right({FromFrame::origin, Velocity<FromFrame>()})).velocity());
}

}  // namespace physics
}  // namespace principia
