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
  // same size.  These functions are fairly expensive.
  std::map<Time, Vector<Length, Frame>> Positions() const;
  std::map<Time, Vector<Speed, Frame>> Velocities() const;
  std::list<Time> Times() const;

  // Return the most recent position/velocity/time.  These functions are dirt
  // cheap.
  Vector<Length, Frame> const& last_position() const;
  Vector<Speed, Frame> const& last_velocity() const;
  Time const& last_time() const;

  // Appends one point to the trajectory.
  void Append(Vector<Length, Frame> const& position,
              Vector<Speed, Frame> const& velocity,
              Time const& time);

  // Removes all data for times (strictly) greater than |time|, as well as all
  // child trajectories forked at times (strictly) greater than |time|.
  void ForgetAfter(Time const& time);

  // Removes all data for times less than or equal to |time|, as well as all
  // child trajectories forked at times less than or equal to |time|.  This
  // trajectory must be a root.
  void ForgetBefore(Time const& time);

  // Creates a new child trajectory forked at time |time|, and returns it.  The
  // child trajectory may be changed independently from the parent trajectory
  // for any time (strictly) greater than |time|.  The child trajectory is owned
  // by its parent trajectory.  Calling ForgetAfter or ForgetBefore on the
  // parent trajectory with an argument that causes the time |time| to be
  // removed deletes the child trajectory.  Deleting the parent trajectory
  // deletes all child trajectories.  |time| must be one of the times of the
  // current trajectory (as returned by Times()).  No transfer of ownership.
  Trajectory* Fork(Time const& time);

  // Returns true if this is a root trajectory.
  bool is_root() const;

  // Returns the root trajectory.
  Trajectory const* root() const;
  Trajectory* root();

  // Returns the fork time for a nonroot trajectory and null for a root
  // trajectory.
  Time const* fork_time() const;

  // The body to which this trajectory pertains.  No transfer of ownership.
  Body const* body() const;

 private:
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

  // A constructor for creating a child trajectory during forking.
  Trajectory(Body const* const body,
             Trajectory* const parent,
             typename States::iterator const& parent_state);

  template<typename Value>
  std::map<Time, Value> GetState(std::function<Value(State const&)> fun) const;

  Body const* const body_;  // Never null.

  Trajectory* const parent_;  // Null for a root trajectory.

  // Null for a root trajectory.
  std::unique_ptr<typename States::iterator> parent_state_;

  // There may be several forks starting from the same time, hence the multimap.
  // Child trajectories are owned.
  std::multimap<Time, std::unique_ptr<Trajectory>> children_;

  States states_;
};

}  // namespace physics
}  // namespace principia

#include "trajectory_body.hpp"
