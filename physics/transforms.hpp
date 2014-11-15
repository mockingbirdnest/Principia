#pragma once

#include "physics/trajectory.hpp"

namespace principia {
namespace physics {

//TODO(phl):Comment
template<typename FromFrame, typename ThroughFrame, typename ToFrame>
class Transforms {
 public:
  static Transforms BodyCentredNonRotating(
      Trajectory<FromFrame> const& centre_trajectory);
  static Transforms BarycentricRotating(
      Trajectory<FromFrame> const& primary_trajectory,
      Trajectory<FromFrame> const& secondary_trajectory);

  typename Trajectory<FromFrame>::template TransformingIterator<ThroughFrame>
  first(Trajectory<FromFrame> const* from_trajectory);

  typename Trajectory<ThroughFrame>:: template TransformingIterator<ToFrame>
  second(Trajectory<ThroughFrame> const* through_trajectory);

 private:
  typename Trajectory<FromFrame>::template Transform<ThroughFrame> first_;

  typename Trajectory<ThroughFrame>::template Transform<ToFrame> second_;
};

}  // namespace physics
}  // namespace principia

#include "physics/transforms_body.hpp"
