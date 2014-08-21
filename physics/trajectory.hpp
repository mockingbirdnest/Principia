#pragma once

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <set>

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::Vector;
using principia::quantities::Acceleration;
using principia::quantities::Length;
using principia::quantities::Speed;
using principia::quantities::Time;

namespace principia {
namespace physics {

class Body;

template<typename Frame>
class Trajectory {
 public:
  // No transfer of ownership.
  explicit Trajectory(Body const* body);
  ~Trajectory() = default;

  // These functions return the series of positions/velocities/times for the
  // trajectory of the body.  All three containers are guaranteed to have the
  // same size.
  std::map<Time, Vector<Length, Frame>> Positions() const;
  std::map<Time, Vector<Speed, Frame>> Velocities() const;
  std::list<Time> Times() const;

  // Appends one point to the trajectory.
  //TODO(phl):Dirtying?
  void Append(Vector<Length, Frame> const& position,
              Vector<Speed, Frame> const& velocity,
              Time const& time);

  // Removes all data for times (strictly) greater than |time|, as well as all
  // child trajectories forked at times (strictly) greater than |time|.  Bursts
  // that start after |time| are also removed.
  void ForgetAfter(Time const& time);

  // Removes all data for times less than or equal to |time|, as well as all
  // child trajectories forked at times less than or equal to |time|.
  //TODO(phl): Bursts are messy.
  void ForgetBefore(Time const& time);

  // Creates a new child trajectory forked at time |time|, and returns it.  The
  // child trajectory may be changed independently from the parent trajectory
  // for any time (strictly) greater than |time|.  Dirtying the parent
  // trajectory before the time |time| dirties the child trajectory.  Deleting
  // the parent trajectory deletes all child trajectories.  |time| must be one
  // of the times of the current trajectory (as returned by times()).  No
  // transfer of ownership.
  Trajectory* Fork(Time const& time);

  // The body to which this trajectory pertains.  No transfer of ownership.
  Body const* body() const;

  // The acceleration applies over the given interval.  It represents e.g., the
  // intrinsic acceleration of the engine.
  void AddBurst(Vector<Acceleration, Frame> const& acceleration,
                Time const& time1,
                Time const& time2);

 private:
  class Burst {
   public:
    Burst(Vector<Acceleration, Frame> const& acceleration,
          Time const& duration);
   private:
    Vector<Acceleration, Frame> const acceleration_;
    Time const duration_;
  };

  class State {
   public:
    State(Vector<Length, Frame> const& position,
          Vector<Speed, Frame> const& velocity);
    Vector<Length, Frame> const& position() const;
    Vector<Speed, Frame> const& velocity() const;
   private:
    Vector<Length, Frame> const position_;
    Vector<Speed, Frame> const velocity_;
  };

  typedef std::map<Time, State> States;

  Trajectory(Body const* const body,
             Trajectory const* const parent,
             typename States::iterator const& parent_state);

  template<typename Value>
  std::map<Time, Value> GetState(std::function<Value(State const&)> fun) const;

  Body const* const body_;  // Never null.

  Trajectory const* const parent_;  // Null for a root trajectory.

  // Null for a root trajectory.
  std::unique_ptr<typename States::iterator> parent_state_;

  // There may be several forks starting from the same time, hence the multimap.
  // Child trajectories are owned.
  std::multimap<Time, std::unique_ptr<Trajectory>> children_;

  States states_;

  std::map<Time, Burst> bursts_;
};

}  // namespace physics
}  // namespace principia

#include "trajectory_body.hpp"
