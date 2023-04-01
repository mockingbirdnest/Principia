#pragma once

#include "physics/similar_motion.hpp"

namespace principia {
namespace physics {
namespace _similar_motion {
namespace internal {

template<typename FromFrame, typename ToFrame>
SimilarMotion<FromFrame, ToFrame>::SimilarMotion(
    Similarity<FromFrame, ToFrame> const& similarity,
    AngularVelocity<FromFrame> const& angular_velocity_of_to_frame,
    Velocity<FromFrame> const& velocity_of_to_frame_origin,
    Variation<double> const& dilatation_rate)
    : similarity_(similarity),
      angular_velocity_of_to_frame_(angular_velocity_of_to_frame),
      velocity_of_to_frame_origin_(velocity_of_to_frame_origin),
      dilatation_rate_(dilatation_rate) {}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
SimilarMotion<FromFrame, ToFrame>::SimilarMotion(
    Similarity<FromFrame, ToFrame> const& similarity,
    AngularVelocity<FromFrame> const& angular_velocity_of_from_frame,
    Velocity<FromFrame> const& velocity_of_from_frame_origin,
    Variation<double> const& dilatation_rate)
    : similarity_(similarity),
      angular_velocity_of_to_frame_(angular_velocity_of_to_frame),
      velocity_of_to_frame_origin_(velocity_of_to_frame_origin),
      dilatation_rate_(dilatation_rate) {}


template<typename FromFrame, typename ToFrame>
DegreesOfFreedom<ToFrame> SimilarMotion<FromFrame, ToFrame>::operator()(
    DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const {
  return {similarity_(degrees_of_freedom.position()),
          orthogonal_map()(
              degrees_of_freedom.velocity() - velocity_of_to_frame_origin_ -
              angular_velocity_of_to_frame_ *
                  (degrees_of_freedom.position() -
                   rigid_transformation_.Inverse()(ToFrame::origin)) /
                  Radian)};
}

}  // namespace internal
}  // namespace _similar_motion
}  // namespace physics
}  // namespace principia
