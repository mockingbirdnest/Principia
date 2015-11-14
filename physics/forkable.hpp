#pragma once

#include <deque>
#include <experimental/optional>  // NOLINT
#include <map>
#include <memory>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {

using geometry::Instant;

namespace physics {

// Forkable and ForkableIterator both use CRTP to achieve static polymorphism on
// the parameters and return type of the member functions: we want them to
// return Tr4jectory and It3rator, not Forkable and ForkableIterator, so that
// the clients don't have to down_cast or construct objects of subclasses.
// ForkableIterator is seen by the clients as a class nested within Forkable.
// However, this cannot be implemented that way because the two classes are
// mutually dependent.  Instead we have two distinct classes: ForkableIterator
// must be instantiated first using an incomplete declaration of Forkable,
// and Forkable may then be instantiated using ForkableIterator.
// The template parameters with 1337 names are those that participate in this
// mutual CRTP.

template<typename Tr4jectory, typename It3rator>
class Forkable;

namespace internal {

// This traits class must export declarations similar to the following:
//
// using TimelineConstIterator = ...;
// static Instant const& time(TimelineConstIterator const it);
//
// TimelineConstIterator must be an STL-like iterator in the timeline of
// Tr4jectory.  |time()| must return the corresponding time.
//
// NOTE(phl): This was originally written as a trait under the assumption that
// we would want to expose STL iterators to clients.  This doesn't seem like a
// good idea anymore, so maybe this should turn into another CRTP class.
template<typename Tr4jectory>
struct ForkableTraits;

// A template for iterating over the timeline of a Forkable object, taking forks
// into account.
template<typename Tr4jectory, typename It3rator>
class ForkableIterator {
  using TimelineConstIterator =
      typename ForkableTraits<Tr4jectory>::TimelineConstIterator;

 public:
  ForkableIterator() = default;

  bool operator==(It3rator const& right) const;
  bool operator!=(It3rator const& right) const;

  It3rator& operator++();
  It3rator& operator--();

 protected:
  // The API that must be implemented by subclasses.
  // Must return |this| of the proper type.
  virtual not_null<It3rator*> that() = 0;
  virtual not_null<It3rator const*> that() const = 0;

  // Returns the point in the timeline that is denoted by this iterator.
  TimelineConstIterator current() const;

 private:
  // Returns the (most forked) trajectory to which this iterator applies.
  not_null<Tr4jectory const*> trajectory() const;

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

  template<typename, typename>
  friend class physics::Forkable;
};

}  // namespace internal

// This template represents a trajectory which is forkable and iterable (using
// a ForkableIterator).
template<typename Tr4jectory, typename It3rator>
class Forkable {
 public:
  // An iterator into the timeline of the trajectory.  Must be STL-like.
  // Beware, if these iterators are invalidated all the guarantees of Forkable
  // are void.
  using TimelineConstIterator =
      typename internal::ForkableTraits<Tr4jectory>::TimelineConstIterator;

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

  It3rator Begin() const;
  It3rator End() const;

  It3rator Find(Instant const& time) const;
  It3rator LowerBound(Instant const& time) const;

  // Returns an iterator denoting the fork point of this object.  Fails if this
  // object is a root.
  It3rator Fork() const;

  // Returns the number of points in this object.  Complexity is O(|length| +
  // |depth|).
  int Size() const;

  void WritePointerToMessage(
      not_null<serialization::Trajectory::Pointer*> const message) const;

  // |trajectory| must be a root.
  static not_null<Tr4jectory*> ReadPointerFromMessage(
      serialization::Trajectory::Pointer const& message,
      not_null<Tr4jectory*> const trajectory);

 protected:
  // The API that must be implemented by subclasses.

  // Must return |this| of the proper type.
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

  // Creates a new child trajectory forked at the given |timeline_it|, and
  // returns it.  The child trajectory shares its data with the current
  // trajectory for times less than or equal to |timeline_it|.  It may be
  // changed independently from the parent trajectory for any time (strictly)
  // greater than |timeline_it|.  The child trajectory is owned by its parent
  // trajectory.  Deleting the parent trajectory deletes all child trajectories.
  // |timeline_it| may be at end if it denotes the fork time of this object.
  not_null<Tr4jectory*> NewFork(TimelineConstIterator const& timeline_it);

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
  // Constructs an Iterator by wrapping the timeline iterator
  // |position_in_ancestor_timeline| which must be an iterator in the timeline
  // of |ancestor|.  |ancestor| must be an ancestor of this trajectory
  // (it may be this object).  |position_in_ancestor_timeline| may only be at
  // end if it is an iterator in this object (and |ancestor| is this object).
  It3rator Wrap(
      not_null<const Tr4jectory*> const ancestor,
      TimelineConstIterator const position_in_ancestor_timeline) const;

  // There may be several forks starting from the same time, hence the multimap.
  // A level of indirection is needed to avoid referencing an incomplete type in
  // CRTP.
  using Children = std::multimap<Instant, std::unique_ptr<Tr4jectory>>;

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

  template<typename, typename>
  friend class internal::ForkableIterator;
};

}  // namespace physics
}  // namespace principia

#include "physics/forkable_body.hpp"
