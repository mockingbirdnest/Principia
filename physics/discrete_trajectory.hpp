#pragma once

#include <functional>
#include <list>
#include <map>
#include <memory>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/forkable.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {

using base::not_null;
using geometry::Instant;
using geometry::Vector;
using geometry::Velocity;
using quantities::Acceleration;
using quantities::Length;
using quantities::Speed;

namespace physics {

template<typename Frame>
class DiscreteTrajectory;

template<typename Frame>
struct ForkableTraits<DiscreteTrajectory<Frame>> {
  using TimelineConstIterator =
      typename std::map<Instant, DegreesOfFreedom<Frame>>::const_iterator;
  static Instant const& time(TimelineConstIterator const it);
};

template<typename Frame>
class DiscreteTrajectory : public Forkable<DiscreteTrajectory<Frame>> {
  using Timeline = std::map<Instant, DegreesOfFreedom<Frame>>;
  using TimelineConstIterator =
      typename Forkable<DiscreteTrajectory<Frame>>::TimelineConstIterator;

 public:
  class Iterator;

  DiscreteTrajectory() = default;
  ~DiscreteTrajectory() override;

  DiscreteTrajectory(DiscreteTrajectory const&) = delete;
  DiscreteTrajectory(DiscreteTrajectory&&) = delete;
  DiscreteTrajectory& operator=(DiscreteTrajectory const&) = delete;
  DiscreteTrajectory& operator=(DiscreteTrajectory&&) = delete;

  // Sets a callback to be run before this trajectory gets destroyed.
  void set_on_destroy(
      std::function<void(not_null<DiscreteTrajectory<Frame>const*> const)>
          on_destroy);

  // TODO(phl): Many/most of the iterator functions are obsolete.  Remove them
  // and use the ones from Forkable.

  // Returns an iterator at the last point of the trajectory.  Complexity is
  // O(1).  The trajectory must not be empty.
  Iterator last() const;

  // These functions return the series of positions/velocities/times for the
  // trajectory.  All three containers are guaranteed to have the same size.
  // These functions are O(|depth| + |length|).
  std::map<Instant, Position<Frame>> Positions() const;
  std::map<Instant, Velocity<Frame>> Velocities() const;
  std::list<Instant> Times() const;

  // Creates a new child trajectory forked at time |time|, and returns it.  The
  // child trajectory shares its data with the current trajectory for times less
  // than or equal to |time|, and is an exact copy of the current trajectory for
  // times greater than |time|.  It may be changed independently from the
  // parent trajectory for any time (strictly) greater than |time|.  The child
  // trajectory is owned by its parent trajectory.  Deleting the parent
  // trajectory deletes all child trajectories.  |time| must be one of the times
  // of this trajectory, and must be at or after the fork time, if any.
  not_null<DiscreteTrajectory<Frame>*> NewForkWithCopy(Instant const& time);

  // Appends one point to the trajectory.
  void Append(Instant const& time,
              DegreesOfFreedom<Frame> const& degrees_of_freedom);

  // Removes all data for times (strictly) greater than |time|, as well as all
  // child trajectories forked at times (strictly) greater than |time|.  |time|
  // must be at or after the fork time, if any.
  void ForgetAfter(Instant const& time);

  // Removes all data for times less than or equal to |time|, as well as all
  // child trajectories forked at times less than or equal to |time|.  This
  // trajectory must be a root.
  void ForgetBefore(Instant const& time);

  // This trajectory must be a root.
  void WriteToMessage(not_null<serialization::Trajectory*> const message) const;

  // NOTE(egg): This should return a |not_null|, but we can't do that until
  // |not_null<std::unique_ptr<T>>| is convertible to |std::unique_ptr<T>|, and
  // that requires a VS 2015 feature (rvalue references for |*this|).
  static std::unique_ptr<DiscreteTrajectory> ReadFromMessage(
      serialization::Trajectory const& message);

  class Iterator : public Forkable<DiscreteTrajectory<Frame>>::Iterator {
   public:
    Iterator(typename Forkable<DiscreteTrajectory<Frame>>::Iterator it);

    Instant const& time() const;
    DegreesOfFreedom<Frame> const& degrees_of_freedom() const;

   private:
    friend class DiscreteTrajectory;
  };

 protected:
  // The API inherited from Forkable.
  not_null<DiscreteTrajectory*> that() override;
  not_null<DiscreteTrajectory const*> that() const override;

  TimelineConstIterator timeline_begin() const override;
  TimelineConstIterator timeline_end() const override;
  TimelineConstIterator timeline_find(Instant const& time) const override;
  TimelineConstIterator timeline_lower_bound(
                            Instant const& time) const override;
  bool timeline_empty() const override;

 private:
  // This trajectory need not be a root.
  void WriteSubTreeToMessage(
      not_null<serialization::Trajectory*> const message) const;

  void FillSubTreeFromMessage(serialization::Trajectory const& message);

  Timeline timeline_;

  std::function<void(not_null<DiscreteTrajectory<Frame>const *> const)>
      on_destroy_;

  template<typename Tr4jectory>
  friend class Forkable;

  // For using the private constructor in maps.
  template<typename, typename>
  friend struct std::pair;
};

}  // namespace physics
}  // namespace principia

#include "discrete_trajectory_body.hpp"
