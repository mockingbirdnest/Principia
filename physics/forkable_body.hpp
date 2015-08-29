#include "forkable.hpp"
#pragma once

#include "physics/forkable.hpp"

namespace principia {
namespace physics {

template<typename Tr4jectory, typename TimelineConstIterator_>
Forkable<Tr4jectory, TimelineConstIterator_>::Forkable() : parent_(nullptr) {}
//TODO(phl): And the other fields?

template<typename Tr4jectory, typename TimelineConstIterator_>
not_null<Tr4jectory*>
Forkable<Tr4jectory, TimelineConstIterator_>::NewFork(Instant const & time) {
  CHECK(timeline_find(time) != timeline_end() ||
        (!is_root() && time == position_in_parent_children_->first))
      << "NewFork at nonexistent time " << time;

  // May be at |timeline_end()|.
  auto timeline_it = timeline_find(time);

  // First create a child in the multimap.
  auto const child_it = children_.emplace(time, Tr4jectory());

  // Now set the members of the child object.
  auto& child_forkable = child_it->second;
  child_forkable.parent_ = that();
  child_forkable.position_in_parent_children_ = child_it;
  child_forkable.position_in_parent_timeline_ = timeline_it;

  // Copy the tail of the trajectory in the child object.
  if (timeline_it != timeline_end()) {
    child_forkable.timeline_insert(++timeline_it, timeline_end());
  }
  return &child_forkable;
}

template<typename Tr4jectory, typename TimelineConstIterator_>
void Forkable<Tr4jectory, TimelineConstIterator_>::DeleteFork(
    not_null<Forkable**> const forkable) {
  CHECK_NOTNULL(*forkable);
  Instant const* const fork_time = (*forkable)->ForkTime();
  CHECK_NOTNULL(fork_time);
  // Find the position of |*forkable| among our children and remove it.
  auto const range = children_.equal_range(*fork_time);
  for (auto it = range.first; it != range.second; ++it) {
    if (&it->second == *forkable) {
      children_.erase(it);
      *forkable = nullptr;
      return;
    }
  }
  LOG(FATAL) << "argument is not a child of this forkable";
}

template<typename Tr4jectory, typename TimelineConstIterator_>
bool Forkable<Tr4jectory, TimelineConstIterator_>::is_root() const {
  return parent_ == nullptr;
}

template<typename Tr4jectory, typename TimelineConstIterator_>
not_null<Tr4jectory const*>
Forkable<Tr4jectory, TimelineConstIterator_>::root() const {
  Tr4jectory const* ancestor = that();
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Tr4jectory, typename TimelineConstIterator_>
not_null<Tr4jectory*>
Forkable<Tr4jectory, TimelineConstIterator_>::root() {
  Tr4jectory* ancestor = that();
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Tr4jectory, typename TimelineConstIterator_>
Instant const* Forkable<Tr4jectory, TimelineConstIterator_>::ForkTime() const {
  if (is_root()) {
    return nullptr;
  } else {
    return &position_in_parent_children_->first;
  }
}

template<typename Tr4jectory, typename TimelineConstIterator_>
bool Forkable<Tr4jectory, TimelineConstIterator_>::Iterator::operator==(
    Iterator const & right) const {
  bool const this_at_end = at_end();
  bool const right_at_end = right.at_end();
  if (this_at_end != right_at_end) {
    return false;
  } else if (this_at_end) {
    return true;
  } else {
    return ancestry_ == right.ancestry_ && current_ == right.current_;
  }
}

template<typename Tr4jectory, typename TimelineConstIterator_>
bool Forkable<Tr4jectory, TimelineConstIterator_>::Iterator::operator!=(
    Iterator const & right) const {
  return !(*this == right);
}

template<typename Tr4jectory, typename TimelineConstIterator_>
typename Forkable<Tr4jectory, TimelineConstIterator_>::Iterator&
Forkable<Tr4jectory, TimelineConstIterator_>::Iterator::operator++() {
  CHECK(!ancestry_.empty());
  CHECK(!at_end());

  Instant const& current_time = *current_/*->first*/;///Traits?

  // Check if there is a next child in the ancestry.
  auto ancestry_it = ancestry_.begin();
  if (++ancestry_it != ancestry_.end()) {
    // There is a next child.  See if we reached its fork time.
    Instant const& current_time = *current_/*->first*/;///Traits?
    not_null<Tr4jectory const*> child = *ancestry_it;
    Instant const& child_fork_time = child->position_in_parent_children_->first;
    if (current_time == child_fork_time) {
      // We have reached the fork time of the next child.  Drop the leading
      // ancestor.
      ancestry_.pop_front();
      // We'd like to iterate over the timeline of the next child, but that
      // timeline may be empty.  So we must skip any empty timeline until either
      // we find a non-empty one or we reach the last child.  All children with
      // empty timelines should fork at the same time.
      while (child->timeline_empty() && ++ancestry_it != ancestry_.end()) {
        child = *ancestry_it;
        ancestry_.pop_front();
        CHECK_EQ(child_fork_time, child->position_in_parent_children_->first);
      }
      // Start iterating over the next child timeline.
      current_ = child->timeline_begin();
      return *this;
    }
  }
  // Business as usual, keep moving along the same timeline.
  ++current_;

  return *this;
}

template<typename Tr4jectory, typename TimelineConstIterator_>
typename Forkable<Tr4jectory, TimelineConstIterator_>::TimelineConstIterator
Forkable<Tr4jectory, TimelineConstIterator_>::Iterator::current() const {
  return current_;
}

template<typename Tr4jectory, typename TimelineConstIterator_>
bool Forkable<Tr4jectory, TimelineConstIterator_>::Iterator::at_end() const {
  return current_ == ancestry_.front()->timeline_end();
}

template<typename Tr4jectory, typename TimelineConstIterator_>
typename Forkable<Tr4jectory, TimelineConstIterator_>::Iterator
Forkable<Tr4jectory, TimelineConstIterator_>::End() const {
  Iterator iterator;
  not_null<Tr4jectory const*> ancestor = that();
  iterator.ancestry_.push_front(ancestor);
  iterator.current_ = ancestor->timeline_end();
  return iterator;
}

template<typename Tr4jectory, typename TimelineConstIterator_>
typename Forkable<Tr4jectory, TimelineConstIterator_>::Iterator
Forkable<Tr4jectory, TimelineConstIterator_>::Find(Instant const& time) const {
  Iterator iterator;

  // Go up the ancestry chain until we find a timeline that covers |time| (that
  // is, |time| is after the first time of the timeline).  Set |current_| to
  // the location of |time|, which may be |end()|.  The ancestry has |forkable|
  // at the back, and the object containing |current_| at the front.
  not_null<Tr4jectory const*> ancestor = that();
  do {
    iterator.ancestry_.push_front(ancestor);
    if (!ancestor->timeline_empty() &&
        *ancestor->timeline_begin()/*->first*/ <= time) {
      iterator.current_ = ancestor->timeline_find(time);  // May be at end.
      break;
    }
    iterator.current_ = ancestor->timeline_end();
    ancestor = ancestor->parent_;
  } while (ancestor != nullptr);

  return iterator;
}


template<typename Tr4jectory, typename TimelineConstIterator_>
typename Forkable<Tr4jectory, TimelineConstIterator_>::Iterator
Forkable<Tr4jectory, TimelineConstIterator_>::Wrap(
    not_null<const Tr4jectory*> const ancestor,
    TimelineConstIterator const position_in_ancestor_timeline) const {
  Iterator iterator;

  // Go up the ancestry chain until we find |ancestor| and set |current_| to
  // |position_in_ancestor_timeline|.  The ancestry has |forkable|
  // at the back, and the object containing |current_| at the front.  If
  // |ancestor| is not found in the ancestry, stop at the root.
  not_null<Tr4jectory const*> ancest0r = that();
  do {
    iterator.ancestry_.push_front(ancest0r);
    if (ancestor == ancest0r) {
      iterator.current_ = position_in_ancestor_timeline;  // May be at end.
      break;
    }
    iterator.current_ = ancest0r->timeline_end();
    ancest0r = ancest0r->parent_;
  } while (ancest0r != nullptr);
  return iterator;
}

}  // namespace physics
}  // namespace principia
