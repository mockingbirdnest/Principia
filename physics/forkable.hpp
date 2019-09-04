
#pragma once

#include <deque>
#include <optional>
#include <map>
#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_forkable {

using base::not_null;
using geometry::Instant;

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

// This traits class must export declarations similar to the following:
//
// using TimelineConstIterator = ...;
// static Instant const& time(TimelineConstIterator it);
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
  virtual ~ForkableIterator() = default;

  // Returns the (most forked) trajectory to which this iterator applies.
  not_null<Tr4jectory const*> trajectory() const;

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
  TimelineConstIterator const& current() const;

 private:
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
  friend class Forkable;
};

// This template represents a trajectory which is forkable and iterable (using
// a ForkableIterator).
template<typename Tr4jectory, typename It3rator>
class Forkable {
 public:
  // An iterator into the timeline of the trajectory.  Must be STL-like.
  // Beware, if these iterators are invalidated all the guarantees of Forkable
  // are void.
  using TimelineConstIterator =
      typename ForkableTraits<Tr4jectory>::TimelineConstIterator;

  Forkable() = default;
  virtual ~Forkable() = default;

  // Cannot be moved or copied because of parent/children pointers.
  Forkable(Forkable const&) = delete;
  Forkable(Forkable&&) = delete;
  Forkable& operator=(Forkable const&) = delete;
  Forkable& operator=(Forkable&&) = delete;

  // Deletes the child trajectory denoted by |trajectory|, which must be a
  // pointer previously returned by NewFork for this object.  Nulls
  // |trajectory|.
  void DeleteFork(Tr4jectory*& trajectory);

  // Returns true if this is a root trajectory.
  bool is_root() const;

  // Returns the root trajectory.
  not_null<Tr4jectory const*> root() const;
  not_null<Tr4jectory*> root();

  not_null<Tr4jectory const*> parent() const;
  not_null<Tr4jectory*> parent();

  It3rator begin() const;
  It3rator end() const;

  typename It3rator::reference front() const;
  typename It3rator::reference back() const;

  It3rator Find(Instant const& time) const;
  It3rator LowerBound(Instant const& time) const;

  // Returns an iterator denoting the fork point of this object.  Fails if this
  // object is a root.
  It3rator Fork() const;

  // Returns the number of points in this object.  Complexity is O(|length| +
  // |depth|).
  std::int64_t Size() const;

  // Returns true if this object is empty.  Complexity is O(1).
  bool Empty() const;

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
  virtual std::int64_t timeline_size() const = 0;

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

  // |fork| must be a non-empty root and its first point must be at the same
  // time as the last point of this object.  |fork| is attached to this object
  // as a child at the end of the timeline.  The caller must then delete the
  // first point of |fork|'s timeline.
  void AttachForkToCopiedBegin(not_null<std::unique_ptr<Tr4jectory>> fork);

  // This object must not be a root.  It is detached from its parent and becomes
  // a root.  All the children which were fork at this object's fork time are
  // changed to be forked at the beginning of this object's timeline.  This
  // requires the caller to ensure that this object's timeline is not empty and
  // that its beginning properly represents the fork time.  Returns an owning
  // pointer to this object.
  not_null<std::unique_ptr<Tr4jectory>> DetachForkWithCopiedBegin();

  // Deletes all forks for times (strictly) greater than |time|.  |time| must be
  // at or after the fork time of this trajectory, if any.
  void DeleteAllForksAfter(Instant const& time);

  // Checks that there exist no forks for times (strictly) less than |time|.
  // This trajectory must be a root.
  void CheckNoForksBefore(Instant const& time);

  // This trajectory need not be a root.  As forks are encountered during tree
  // traversal their pointer is nulled-out in |forks|.
  void WriteSubTreeToMessage(
      not_null<serialization::DiscreteTrajectory*> message,
      std::vector<Tr4jectory*>& forks) const;

  void FillSubTreeFromMessage(serialization::DiscreteTrajectory const& message,
                              std::vector<Tr4jectory**> const& forks);

 private:
  // Constructs an Iterator by wrapping the timeline iterator
  // |position_in_ancestor_timeline| which must be an iterator in the timeline
  // of |ancestor|.  |ancestor| must be an ancestor of this trajectory
  // (it may be this object).  |position_in_ancestor_timeline| may only be at
  // end if it is an iterator in this object (and |ancestor| is this object).
  It3rator Wrap(not_null<Tr4jectory const*> ancestor,
                TimelineConstIterator position_in_ancestor_timeline) const;

  // There may be several forks starting from the same time, hence the multimap.
  // A level of indirection is needed to avoid referencing an incomplete type in
  // CRTP.
  using Children = std::multimap<Instant, std::unique_ptr<Tr4jectory>>;

  // Null for a root.
  Tr4jectory* parent_ = nullptr;

  // This iterator is never at |end()|.
  std::optional<typename Children::iterator> position_in_parent_children_;

  // This iterator is at |end()| if the fork time is not in the parent timeline,
  // i.e. is the parent timeline's own fork time.
  std::optional<TimelineConstIterator> position_in_parent_timeline_;
  Children children_;

  template<typename, typename>
  friend class ForkableIterator;
};

}  // namespace internal_forkable

using internal_forkable::Forkable;

}  // namespace physics
}  // namespace principia

#include "physics/forkable_body.hpp"
