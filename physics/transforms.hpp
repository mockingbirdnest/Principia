#pragma once

#include "physics/trajectory.hpp"

namespace principia {
namespace physics {

template<typename FromFrame, typename ToFrame>
Trajectory<FromFrame>::TransformingIterator<ToFrame>
BodyCentredNonRotatingTransformingIterator(
    Trajectory<FromFrame> const& centre_trajectory,
    Trajectory<FromFrame> const* transformed_trajectory);

}  // namespace physics
}  // namespace principia

#include "physics/transforms_body.hpp"
