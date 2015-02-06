#pragma once

#include <functional>
#include <list>
#include <map>
#include <memory>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/physics.pb.h"

using principia::base::not_null;
using principia::geometry::Instant;
using principia::geometry::Vector;
using principia::geometry::Velocity;
using principia::quantities::Acceleration;
using principia::quantities::Length;
using principia::quantities::Speed;

namespace principia {
namespace physics {

class Body;

template<typename Frame>
class Trajectory {
  // There may be several forks starting from the same time, hence the multimap.
  using Children =
      std::multimap<Instant, not_null<std::unique_ptr<Trajectory>>>;
  using Timeline = std::map<Instant, DegreesOfFreedom<Frame>>;

  // The two iterators denote entries in the containers of the parent, and they
  // are never past the end.  Therefore, they are not invalidated by swapping
  // the containers of the parent.
  struct Fork {
    typename Children::const_iterator children;
    typename Timeline::const_iterator timeline;
  };

 public:
  class NativeIterator;
  template<typename ToFrame>
  class TransformingIterator;

  // A function that transforms the coordinates to a different frame.
  template<typename ToFrame>
  using Transform = std::function<DegreesOfFreedom<ToFrame>(
                        Instant const&,
                        DegreesOfFreedom<Frame> const&,
                        not_null<Trajectory<Frame> const*> const)>;

  // No transfer of ownership.  |body| must live longer than the trajectory as
  // the trajectory holds a reference to it.  If |body| is oblate it must be
  // expressed in the same frame as the trajectory.
  explicit Trajectory(not_null<Body const*> const body);
  ~Trajectory() = default;

  Trajectory(Trajectory const&) = delete;
  Trajectory(Trajectory&&) = delete;
  Trajectory& operator=(Trajectory const&) = delete;
  Trajectory& operator=(Trajectory&&) = delete;

  // Returns an iterator at the first point of the trajectory.  Complexity is
  // O(|depth|).  The result may be at end if the trajectory is empty.
  NativeIterator first() const;

  // Returns at the first point of the trajectory which is on or after |time|.
  // Complexity is O(|depth| + Ln(|length|)).  The result may be at end if the
  // |time| is after the end of the trajectory.
  NativeIterator on_or_after(Instant const& time) const;

  // Returns an iterator at the last point of the trajectory.  Complexity is
  // O(1).  The trajectory must not be empty.
  NativeIterator last() const;

  // Same as |first| above, but returns an iterator that performs a coordinate
  // tranformation to ToFrame.
  template<typename ToFrame>
  TransformingIterator<ToFrame> first_with_transform(
      Transform<ToFrame> const& transform) const;

  // Returns at the first point of the trajectory which is on or after |time|.
  // Complexity is O(|depth| + Ln(|length|)).  The result may be at end if the
  // |time| is after the end of the trajectory.
  template<typename ToFrame>
  TransformingIterator<ToFrame> on_or_after_with_transform(
      Instant const& time,
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
  // child trajectory shares its data with the current trajectory for times less
  // than or equal to |time|, and is an exact copy of the current trajectory for
  // times greater than |time|.  It may be changed independently from the
  // parent trajectory for any time (strictly) greater than |time|.  The child
  // trajectory is owned by its parent trajectory.  Calling ForgetAfter or
  // ForgetBefore on the parent trajectory with an argument that causes the time
  // |time| to be removed deletes the child trajectory.  Deleting the parent
  // trajectory deletes all child trajectories.  |time| must be one of the times
  // of the current trajectory (as returned by Times()).  No transfer of
  // ownership.
  not_null<Trajectory*> NewFork(Instant const& time);

  // Deletes the child trajectory denoted by |*fork|, which must be a pointer
  // previously returned by NewFork for this object.  Nulls |*fork|.
  void DeleteFork(not_null<Trajectory**> const fork);

  // Returns true if this is a root trajectory.
  bool is_root() const;

  // Returns the root trajectory.
  not_null<Trajectory const*> root() const;
  not_null<Trajectory*> root();

  // Returns the fork time for a nonroot trajectory and null for a root
  // trajectory.
  Instant const* fork_time() const;

  // The body to which this trajectory pertains.  The body is cast to the type
  // B.  An error occurs in debug mode if the cast fails.
  template<typename B>
  std::enable_if_t<std::is_base_of<Body, B>::value,
                   not_null<B const*>> body() const;

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

  // This trajectory must be a root.  The intrinsic acceleration is not
  // serialized.  The body is not owned, and therefore is not serialized.
  void WriteToMessage(not_null<serialization::Trajectory*> const message) const;

  // NOTE(egg): This should return a |not_null|, but we can't do that until
  // |not_null<std::unique_ptr<T>>| is convertible to |std::unique_ptr<T>|, and
  // that requires a VS 2015 feature (rvalue references for |*this|).
  static std::unique_ptr<Trajectory> ReadFromMessage(
      serialization::Trajectory const& message,
      not_null<Body const*> const body);

  void WritePointerToMessage(
      not_null<serialization::Trajectory::Pointer*> const message) const;

  // |trajectory| must be a root.
  static not_null<Trajectory*> ReadPointerFromMessage(
      serialization::Trajectory::Pointer const& message,
      not_null<Trajectory*> const trajectory);

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
    void InitializeFirst(not_null<Trajectory const*> const trajectory);
    void InitializeOnOrAfter(Instant const& time,
                             not_null<Trajectory const*> const trajectory);
    void InitializeLast(not_null<Trajectory const*> const trajectory);
    typename Timeline::const_iterator current() const;
    not_null<Trajectory const*> trajectory() const;

   private:
    // |ancestry_| has one more element than |forks_|.  The first element in
    // |ancestry_| is the root.  There is no element in |forks_| for the root.
    // It is therefore empty for a root trajectory.
    typename Timeline::const_iterator current_;
    std::list<not_null<Trajectory const*>> ancestry_;  // Pointers not owned.
    std::list<Fork> forks_;
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
    DegreesOfFreedom<ToFrame> degrees_of_freedom() const;
   private:
    explicit TransformingIterator(Transform<ToFrame> const& transform);
    Transform<ToFrame> transform_;
    friend class Trajectory;
  };

 private:
  // A constructor for creating a child trajectory during forking.
  Trajectory(not_null<Body const*> const body,
             not_null<Trajectory*> const parent,
             Fork const& fork);

  // This trajectory need not be a root.
  void WriteSubTreeToMessage(
      not_null<serialization::Trajectory*> const message) const;

  void FillSubTreeFromMessage(serialization::Trajectory const& message);

  not_null<Body const*> const body_;

  // Both of these members are null for a root trajectory.
  std::unique_ptr<Fork> fork_;
  Trajectory* const parent_;

  Children children_;
  Timeline timeline_;

  std::unique_ptr<IntrinsicAcceleration> intrinsic_acceleration_;
};

}  // namespace physics
}  // namespace principia

#include "trajectory_body.hpp"
