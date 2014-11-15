#pragma once

#include "physics/trajectory.hpp"

namespace principia {
namespace physics {

// This class represent a pair of transformations of a trajectory from
// |FromFrame| to |ToFrame| with an intermediate representation in
// |ThroughFrame|.  Both |FromFrame| and |ToFrame| must be inertial frames.
template<typename FromFrame, typename ThroughFrame, typename ToFrame>
class Transforms {
 public:
  // A factory method where |ThroughFrame| is defined as follows: it has the
  // same axes as |FromFrame| and the body of |centre_trajectory| is the origin
  // of |ThroughFrame|.
  static Transforms BodyCentredNonRotating(
      Trajectory<FromFrame> const& from_centre_trajectory,
      Trajectory<ToFrame> const& to_centre_trajectory);

  // A factory method where |ThroughFrame| is defined as follows: its X axis
  // goes from the primary to the secondary bodies, its Y axis is in the plane
  // of the velocities of the bodies in their barycentric frame, on the same
  // side of the X axis as the velocity of the primary body, its Z axis is such
  // that it is right-handed.  The barycentre of the bodies is the origin of
  // |ThroughFrame|.
  static Transforms BarycentricRotating(
      Trajectory<FromFrame> const& from_primary_trajectory,
      Trajectory<ToFrame> const& to_primary_trajectory,
      Trajectory<FromFrame> const& from_secondary_trajectory,
      Trajectory<ToFrame> const& to_secondary_trajectory);

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
