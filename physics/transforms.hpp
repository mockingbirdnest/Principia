#pragma once

#include "physics/trajectory.hpp"

namespace principia {
namespace physics {

// This class represent a pair of transformations of a trajectory from
// |FromFrame| to |ToFrame| with an intermediate representation in
// |ThroughFrame|.
template<typename FromFrame, typename ThroughFrame, typename ToFrame>
class Transforms {
 public:
  // A factory method where the intermediate frame is in translation with the
  // body of |centre_trajectory|.
  static Transforms BodyCentredNonRotating(
      Trajectory<FromFrame> const& centre_trajectory);

  // A factory method where the intermediate frame is in translation with the
  // barycentre of the two bodies, whose X axis goes from the primary to the
  // secondary, whose Y axis is in the plane of the velocities of the bodies in
  // their barycentric frame, along the velocity of the primary body, and which
  // is right-handed.
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
