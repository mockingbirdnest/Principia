#pragma once

#include "physics/trajectory.hpp"

namespace principia {
namespace physics {

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
class Transforms {
 public:
  Transforms BodyCentredNonRotating(
      Trajectory<FromFrame> const& centre_trajectory);

  typename Trajectory<FromFrame>::template TransformingIterator<ThroughFrame>
  first(Trajectory<FromFrame> const* from_trajectory);

  typename Trajectory<ThroughFrame>:: template TransformingIterator<ToFrame>
  second(Trajectory<ThroughFrame> const* through_trajectory);

 private:
  typename Trajectory<FromFrame>::template Transform<ThroughFrame>
  first_transform_;

  typename Trajectory<ThroughFrame>::template Transform<ToFrame>
  second_transform_;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Transforms<FromFrame, ThroughFrame, ToFrame>
BodyCentredNonRotatingTransformingIterator(
    Trajectory<FromFrame> const& centre_trajectory);

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Transforms<FromFrame, ThroughFrame, ToFrame>
BarycentricRotatingTransformingIterator(
    Trajectory<FromFrame> const& primary_trajectory,
    Trajectory<FromFrame> const& secondary_trajectory);

}  // namespace physics
}  // namespace principia

#include "physics/transforms_body.hpp"
