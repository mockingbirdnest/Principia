#pragma once

#include <optional.hpp>

#include <deque>
#include <map>
#include <memory>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {

using geometry::Instant;

namespace physics {

// This traits class must export declarations similar to the following:
//
// using TimelineConstIterator = ...;
// static Instant const& time(TimelineConstIterator const it);
//
// TimelineConstIterator must be an STL-like iterator in the timeline of
// Tr4jectory.  |time()| must return the corresponding time.
template<typename Tr4jectory>
struct ForkableTraits;

// This template represents a trajectory which is forkable and iterable.  It
// uses CRTP to achieve static polymorphism on the return type of the member
// functions: we want them to return Tr4jectory, not Forkable, so that the
// clients don't have to down_cast.
template<typename Tr4jectory>
class Forkable {
 public:
  // An iterator into the timeline of the trajectory.  Must be STL-like.
  // Beware, if these iterators are invalidated all the guarantees of Forkable
  // are void.
  using TimelineConstIterator =
      typename ForkableTraits<Tr4jectory>::TimelineConstIterator;

  Forkable() = default;
  virtual ~Forkable() = default;

  // Deletes the child trajectory denoted by |*trajectory|, which must be a
  // pointer previously returned by NewFork for this object.  Nulls
  // |*trajectory|.
  void DeleteFork(not_null<Tr4jectory**> const trajectory);

  // Returns true if this is a root trajectory.
  bool is_root() const;

  // Returns the root trajectory.
  not_null<Tr4jectory const*> root() const;
  not_null<Tr4jectory*> root();

  // Returns the fork time for a nonroot trajectory and null for a root
  // trajectory.
  std::experimental::optional<Instant> ForkTime() const;

  // A base class for iterating over the timeline of a trajectory, taking forks
  // into account.
  class Iterator {
   public:
    bool operator==(Iterator const& right) const;
    bool operator!=(Iterator const& right) const;

    Iterator& operator++();
    Iterator& operator--();

    // Returns the point in the timeline that is denoted by this iterator.
    TimelineConstIterator current() const;

    // Returns the (most forked) trajectory to which this iterator applies.
    not_null<Tr4jectory const*> trajectory() const;

   private:
    Iterator() = default;

    // We want a single representation for an end iterator.  In various places
    // we may end up with |current_| at the end of its timeline, but that
    // timeline is not the "most forked" one.  This function normalizes this
    // object so that there is only one entry in the |ancestry_| (the "most
    // forked" one) and |current_| is at its end.
    void NormalizeIfEnd();

    // Checks that this object verifies the invariants enforced by
    // NormalizeIfEnd and dies if it doesn't.
    void CheckNormalizedIfEnd();

    // |ancestry_| is never empty.  |current_| is an iterator in the timeline
    // for |ancestry_.front()|.  |current_| may be at end.
    TimelineConstIterator current_;
    std::deque<not_null<Tr4jectory const*>> ancestry_;  // Pointers not owned.

    template<typename Tr4jectory>
    friend class Forkable;
  };

  Iterator Begin() const;
  Iterator End() const;

  Iterator Find(Instant const& time) const;
  Iterator LowerBound(Instant const& time) const;

  // Constructs an Iterator by wrapping the timeline iterator
  // |position_in_ancestor_timeline| which must be an iterator in the timeline
  // of |ancestor|.  |ancestor| must be an ancestor of this trajectory
  // (it may be this object).  |position_in_ancestor_timeline| may only be at
  // end if it is an iterator in this object (and ancestor is this object).
  Iterator Wrap(
      not_null<const Tr4jectory*> const ancestor,
      TimelineConstIterator const position_in_ancestor_timeline) const;

  void WritePointerToMessage(
      not_null<serialization::Trajectory::Pointer*> const message) const;

  // |trajectory| must be a root.
  static not_null<Tr4jectory*> ReadPointerFromMessage(
      serialization::Trajectory::Pointer const& message,
      not_null<Tr4jectory*> const trajectory);

 protected:
  // The API that must be implemented by subclasses.
  // TODO(phl): Try to reduce this API.  Forkable should probably not modify the
  // timeline.

  // Must return |this| of the proper type
  virtual not_null<Tr4jectory*> that() = 0;
  virtual not_null<Tr4jectory const*> that() const = 0;

  // STL-like operations.
  virtual TimelineConstIterator timeline_begin() const = 0;
  virtual TimelineConstIterator timeline_end() const = 0;
  virtual TimelineConstIterator timeline_find(Instant const& time) const = 0;
  virtual TimelineConstIterator timeline_lower_bound(
                                    Instant const& time) const = 0;
  virtual bool timeline_empty() const = 0;

 protected:
  // The API that subclasses may use to implement their public operations.

  // Creates a new child trajectory forked at time |time|, and returns it.  The
  // child trajectory shares its data with the current trajectory for times less
  // than or equal to |time|.  It may be changed independently from the parent
  // trajectory for any time (strictly) greater than |time|.  The child
  // trajectory is owned by its parent trajectory.  Deleting the parent
  // trajectory deletes all child trajectories.  |time| must be one of the times
  // of this trajectory, and must be at or after the fork time, if any.
  not_null<Tr4jectory*> NewFork(Instant const& time);

  // Deletes all forks for times (strictly) greater than |time|.  |time| must be
  // at or after the fork time of this trajectory, if any.
  void DeleteAllForksAfter(Instant const& time);

  // Deletes all forks for times less than or equal to |time|.  This trajectory
  // must be a root.
  void DeleteAllForksBefore(Instant const& time);

  // This trajectory need not be a root.
  void WriteSubTreeToMessage(
      not_null<serialization::Trajectory*> const message) const;

  void FillSubTreeFromMessage(serialization::Trajectory const& message);

 private:
  // There may be several forks starting from the same time, hence the multimap.
  using Children = std::multimap<Instant, Tr4jectory>;

  // Null for a root.
  Tr4jectory* parent_ = nullptr;

  // This iterator is never at |end()|.
  std::experimental::optional<typename Children::const_iterator>
      position_in_parent_children_;

  // This iterator is at |end()| if the fork time is not in the parent timeline,
  // i.e. is the parent timeline's own fork time.
  std::experimental::optional<TimelineConstIterator>
      position_in_parent_timeline_;

  Children children_;
};

}  // namespace physics
}  // namespace principia

#include "physics/forkable_body.hpp"
