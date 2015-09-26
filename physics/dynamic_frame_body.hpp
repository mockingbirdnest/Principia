#pragma once

#include "physics/dynamic_frame.hpp"

namespace principia {



namespace physics {
template <typename FromFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame>::RigidMotion(
    RigidTransformation<FromFrame, ToFrame> const& rigid_transformation,
    AngularVelocity<ToFrae> const& rotation,
    Velocity<ToFrame> const& translation)
    : rigid_transformation_(rigid_transformation),
      rotation_(rotation),
      translation_(translation) {}

template <typename FromFrame, typename ToFrame>
RigidTransformation<FromFrame, ToFrame> const&
RigidMotion<FromFrame, ToFrame>::rigid_transformation() const {
  return rigid_tranformation_;
}

template <typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> const&
RigidMotion<FromFrame, ToFrame>::orthogonal_map() const {
  rigid_tranformation_.linear_map();
}

template <typename FromFrame, typename ToFrame>
DegreesOfFreedom<ToFrame> RigidMotion<FromFrame, ToFrame>::operator()(
    DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const {
  Position<ToFrame> const result_position = rigid_tranformation_(degrees_of_freedom.position);
  return {result_position,
          orthogonal_map()(degrees_of_freedom.velocity) + translation_ +
              rotation_ * (result_position - FromFrame::origin)};
}

template <typename FromFrame, typename ToFrame>
RigidMotion<ToFrame, FromFrame> RigidMotion<FromFrame, ToFrame>::Inverse()
    const {
  return RigidMotion<ToFrame, FromFrame>(rigid_tranformation_.Inverse(),
                                         -orthogonal_map()(rotation_),
                                         -orthogonal_map()(translation_));
}

#ifdef(EGG_IS_NOT_LAZY)
template <typename FromFrame, typename ThroughFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame> operator*(
    RigidMotion<ThroughFrame, ToFrame> const& left,
    RigidMotion<FromFrame, ThroughFrame> const& right) {
  return RigidMotion<FromFrame, ToFrame>();
}
#endif

}  // namespace physics
}  // namespace principia
