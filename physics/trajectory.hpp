#pragma once

#include <functional>
#include <list>
#include <map>
#include <memory>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::Instant;
using principia::geometry::Vector;
using principia::geometry::Velocity;
using principia::quantities::Acceleration;
using principia::quantities::Length;
using principia::quantities::Speed;

namespace principia {
namespace physics {

template<typename Frame>
class Body;

template<typename Frame>
class Trajectory {
 public:
  class NativeIterator;
  template<typename ToFrame>
  class TransformingIterator;

  // A function that transforms the coordinates to a different frame.
  template<typename ToFrame>
  using Transform =
      std::function<DegreesOfFreedom<ToFrame>(Instant const&,
                                              DegreesOfFreedom<Frame> const&)>;

  // No transfer of ownership.  |body| must live longer than the trajectory as
  // the trajectory holds a reference to it.
  explicit Trajectory(Body<Frame> const& body);
  ~Trajectory() = default;

  // Returns an iterator at the first point of the trajectory.  Complexity is
  // O(|depth|).  The result may be at end if the trajectory is empty.
  NativeIterator first() const;

  // Returns an iterator at the last point of the trajectory.  Complexity is
  // O(1).  The trajectory must not be empty.
  NativeIterator last() const;

  // Same as |first| above, but returns an iterator that performs a coordinate
  // tranformation to ToFrame.
  template<typename ToFrame>
  TransformingIterator<ToFrame> first_with_transform(
      Transform<ToFrame> const& transform) const;

  // Same as |last| above, but returns an iterator that performs a coordinate
  // tranformation to ToFrame.
  template<typename ToFrame>
  TransformingIterator<ToFrame> last_with_transform(
      Transform<ToFrame> const& transform) const;

  // These functions return the series of positions/velocities/times for the
  // trajectory of the body.  All three containers are guaranteed to have the
  // same size.  These functions are O(|depth| + |length|).
  std::map<Instant, Position<Frame>> Positions() const;
  std::map<Instant, Velocity<Frame>> Velocities() const;
  std::list<Instant> Times() const;

  // Appends one point to the trajectory.
  void Append(Instant const& time,
              DegreesOfFreedom<Frame> const& degrees_of_freedom);

  // Removes all data for times (strictly) greater than |time|, as well as all
  // child trajectories forked at times (strictly) greater than |time|.
  void ForgetAfter(Instant const& time);

  // Removes all data for times less than or equal to |time|, as well as all
  // child trajectories forked at times less than or equal to |time|.  This
  // trajectory must be a root.
  void ForgetBefore(Instant const& time);

  // Creates a new child trajectory forked at time |time|, and returns it.  The
  // child trajectory may be changed independently from the parent trajectory
  // for any time (strictly) greater than |time|.  The child trajectory is owned
  // by its parent trajectory.  Calling ForgetAfter or ForgetBefore on the
  // parent trajectory with an argument that causes the time |time| to be
  // removed deletes the child trajectory.  Deleting the parent trajectory
  // deletes all child trajectories.  |time| must be one of the times of the
  // current trajectory (as returned by Times()).  No transfer of ownership.
  Trajectory* Fork(Instant const& time);

  // Deletes the child trajectory denoted by |*fork|, which must be a pointer
  // previously returned by Fork for this object.  Nulls |*fork|.
  void DeleteFork(Trajectory** const fork);

  // Returns true if this is a root trajectory.
  bool is_root() const;

  // Returns the root trajectory.
  Trajectory const* root() const;
  Trajectory* root();

  // Returns the fork time for a nonroot trajectory and null for a root
  // trajectory.
  Instant const* fork_time() const;

  // The body to which this trajectory pertains.
  Body<Frame> const& body() const;

  // This function represents the intrinsic acceleration of a body, irrespective
  // of any external field.  It can be due e.g., to an engine burn.
  using IntrinsicAcceleration =
      std::function<Vector<Acceleration, Frame>(Instant const& time)>;

  // Sets the intrinsic acceleration for the trajectory of a massless body.
  // For a nonroot trajectory the intrinsic acceleration only applies to times
  // (strictly) greater than |fork_time()|.  In other words, the function
  // |acceleration| is never called with times less than or equal to
  // |fork_time()|.  It may, however, be called with times beyond |last_time()|.
  // For a root trajectory the intrinsic acceleration applies to times greater
  // than or equal to the first time of the trajectory.  Again, it may apply
  // beyond |last_time()|.
  // It is an error to call this function for a trajectory that already has an
  // intrinsic acceleration, or for the trajectory of a massive body.
  void set_intrinsic_acceleration(IntrinsicAcceleration const acceleration);

  // Removes any intrinsic acceleration for the trajectory.
  void clear_intrinsic_acceleration();

  // Returns true if this trajectory has an intrinsic acceleration.
  bool has_intrinsic_acceleration() const;

  // Computes the intrinsic acceleration for this trajectory at time |time|.  If
  // |has_intrinsic_acceleration()| return false, or if |time| is before the
  // |fork_time()| (or initial time) of this trajectory, the returned
  // acceleration is zero.
  Vector<Acceleration, Frame> evaluate_intrinsic_acceleration(
      Instant const& time) const;

  // A base class for iterating over the timeline of a trajectory, taking forks
  // into account.  Objects of this class cannot be created.
  class Iterator {
   public:
    Iterator& operator++();
    bool at_end() const;
    Instant const& time() const;

   protected:
    using Timeline = std::map<Instant, DegreesOfFreedom<Frame>>;

    Iterator() = default;
    // No transfer of ownership.
    void InitializeFirst(Trajectory const* trajectory);
    void InitializeLast(Trajectory const* trajectory);
    typename Timeline::const_iterator current() const;

   private:
    typename Timeline::const_iterator current_;
    std::list<Trajectory const*> ancestry_;  // Pointers not owned.
    std::list<typename Timeline::iterator> forks_;
  };

  // An iterator which returns the coordinates in the native frame of the
  // trajectory, i.e., |Frame|.
  class NativeIterator : public Iterator {
   public:
    DegreesOfFreedom<Frame> const& degrees_of_freedom() const;
   private:
    NativeIterator() = default;
    friend class Trajectory;
  };

  // An iterator which returns the coordinates in another frame.
  template<typename ToFrame>
  class TransformingIterator : public Iterator {
   public:
    DegreesOfFreedom<ToFrame> const& degrees_of_freedom() const;
   private:
    explicit TransformingIterator(Transform<ToFrame> const& transform);
    Transform<ToFrame> const transform_;
  };

 private:
  using Timeline = std::map<Instant, DegreesOfFreedom<Frame>>;

  // A constructor for creating a child trajectory during forking.
  Trajectory(Body<Frame> const& body,
             Trajectory* const parent,
             typename Timeline::iterator const& fork);

  Body<Frame> const& body_;

  Trajectory* const parent_;  // Null for a root trajectory.

  // Null for a root trajectory.
  std::unique_ptr<typename Timeline::iterator> fork_;

  // There may be several forks starting from the same time, hence the multimap.
  // Child trajectories are owned.
  std::multimap<Instant, std::unique_ptr<Trajectory>> children_;

  Timeline timeline_;

  std::unique_ptr<IntrinsicAcceleration> intrinsic_acceleration_;
};

}  // namespace physics
}  // namespace principia

#include "trajectory_body.hpp"
