#pragma once

#include "physics/trajectory.hpp"

namespace principia {
namespace physics {

template<typename FromFrame, typename ToFrame>
typename Trajectory<FromFrame>::TransformingIterator<ToFrame>
BodyCentredNonRotatingTransformingIterator(
    Trajectory<FromFrame> const& centre_trajectory,
    Trajectory<FromFrame> const* transformed_trajectory);

template<typename FromFrame, typename ToFrame>
typename Trajectory<FromFrame>::TransformingIterator<ToFrame>
BarycentricRotatingTransformingIterator(
    Trajectory<FromFrame> const& primary_trajectory,
    Trajectory<FromFrame> const& secondary_trajectory,
    Trajectory<FromFrame> const* transformed_trajectory);

}  // namespace physics
}  // namespace principia

#include "physics/transforms_body.hpp"
