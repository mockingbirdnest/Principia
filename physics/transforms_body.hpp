#pragma once

#include "glog/logging.h"
#include "physics/transforms.hpp"

namespace principia {
namespace physics {

template<typename FromFrame, typename ToFrame>
typename Trajectory<FromFrame>::TransformingIterator<ToFrame>
BodyCentredNonRotatingTransformingIterator(
    Trajectory<FromFrame> const& centre_trajectory,
    Trajectory<FromFrame> const* transformed_trajectory) {
  CHECK_NOTNULL(transformed_trajectory);
  Trajectory<FromFrame>::Transform<ToFrame> transform =
      [&centre_trajectory](
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom) ->
    DegreesOfFreedom<ToFrame> {
    DegreesOfFreedom<FromFrame> const& last_centre_degrees_of_freedom =
        centre_trajectory.last().degrees_of_freedom();
    // GetDegreesOfFreedom is Ln(N), but it doesn't matter unless the map gets
    // very big, in which case we'll have cache misses anyway.
    DegreesOfFreedom<FromFrame> const& centre_degrees_of_freedom =
        centre_trajectory.GetDegreesOfFreedom(t);
    return {from_degrees_of_freedom.position -
                centre_degrees_of_freedom.position +
                last_centre_degrees_of_freedom.position,
            from_degrees_of_freedom.velocity -
                centre_degrees_of_freedom.velocity};
  };
  return transformed_trajectory->first_with_transform(transform);
}

}  // namespace physics
}  // namespace principia
