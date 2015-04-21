#pragma once

#include <functional>
#include <map>
#include <memory>
#include <utility>

#include "base/not_null.hpp"
#include "physics/trajectory.hpp"

namespace principia {

using base::not_null;

namespace physics {

// This class represent a pair of transformations of a trajectory from
// |FromFrame| to |ToFrame| with an intermediate representation in
// |ThroughFrame|.  Note that the trajectory in |ToFrame| is not the trajectory
// of a body since its past changes from moment to moment.
template<typename Object,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
class Transforms {
  static_assert(FromFrame::is_inertial && ToFrame::is_inertial,
                "Both FromFrame and ToFrame must be inertial");

 public:
  // The trajectories are evaluated lazily because they may be extended or
  // deallocated/reallocated between the time when the transforms are created
  // and the time when they are applied.  Thus, the lambdas couldn't capture the
  // trajectories by value nor by reference.  Instead, they capture a copy of a
  // function that accesses the trajectories.
  //TODO(phl): Fix comment
  template<typename Frame>
  using LazyTrajectory = Trajectory<Frame> const& (Object::*)() const;

  // A factory method where |ThroughFrame| is defined as follows: it has the
  // same axes as |FromFrame| and the body of |centre_trajectory| is the origin
  // of |ThroughFrame|.
  static not_null<std::unique_ptr<Transforms>> BodyCentredNonRotating(
      Object const& centre,
      LazyTrajectory<ToFrame> const& to_trajectory);

  // A factory method where |ThroughFrame| is defined as follows: its X axis
  // goes from the primary to the secondary bodies, its Y axis is in the plane
  // of the velocities of the bodies in their barycentric frame, on the same
  // side of the X axis as the velocity of the primary body, its Z axis is such
  // that it is right-handed.  The barycentre of the bodies is the origin of
  // |ThroughFrame|.
  static not_null<std::unique_ptr<Transforms>> BarycentricRotating(
      Object const& primary,
      Object const& secondary,
      LazyTrajectory<ToFrame> const& to_trajectory);

  // Use this only for testing!
  static not_null<std::unique_ptr<Transforms>> DummyForTesting();

  typename Trajectory<FromFrame>::template TransformingIterator<ThroughFrame>
  first(Object const& object,
        LazyTrajectory<FromFrame> const& from_trajectory);

  typename Trajectory<FromFrame>::template TransformingIterator<ThroughFrame>
  first_on_or_after(Object const& object,
                    LazyTrajectory<FromFrame> const& from_trajectory,
                    Instant const& time);

  typename Trajectory<ThroughFrame>:: template TransformingIterator<ToFrame>
  second(Object const& object,
         LazyTrajectory<ThroughFrame> const& through_trajectory);

 private:
  template<typename Frame1, typename Frame2>
  using LazyTransform = std::function<DegreesOfFreedom<Frame2>(
                            LazyTrajectory<Frame1> const&,
                            Instant const&,
                            DegreesOfFreedom<Frame1> const&,
                            not_null<Trajectory<Frame1> const*> const)>;

  LazyTransform<FromFrame, ThroughFrame> first_;
  typename Trajectory<ThroughFrame>::template Transform<ToFrame> second_;

  // A simple cache with no eviction, which monitors the hit rate.
  template<typename Frame1, typename Frame2>
  class Cache {
   public:
    bool Lookup(not_null<Trajectory<Frame1> const*> const trajectory,
                Instant const& time,
                not_null<DegreesOfFreedom<Frame2>**> degrees_of_freedom);

    void Insert(not_null<Trajectory<Frame1> const*> const trajectory,
                Instant const& time,
                DegreesOfFreedom<Frame2> const& degrees_of_freedom);

   private:
    std::map<std::pair<not_null<Trajectory<Frame1> const*>, Instant const>,
             DegreesOfFreedom<Frame2>> map_;
    std::map<not_null<Trajectory<Frame1> const*>, std::int64_t>
        number_of_lookups_;
    std::map<not_null<Trajectory<Frame1> const*>, std::int64_t> number_of_hits_;
  };

  // A cache for the result of the |first_| transform.  This cache assumes that
  // the iterator is never called with the same time but different degrees of
  // freedom.
  Cache<FromFrame, ThroughFrame> first_cache_;
};

}  // namespace physics
}  // namespace principia

#include "physics/transforms_body.hpp"
