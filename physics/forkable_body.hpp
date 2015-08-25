#include "forkable.hpp"
#pragma once

#include "physics/forkable.hpp"

namespace principia {
namespace physics {

template<typename Tr4jectory>
not_null<Forkable<Tr4jectory>*> Forkable<Tr4jectory>::NewFork(
    Instant const & time) {
  CHECK(timeline_find(time) != timeline_end() ||
        (!is_root() && time == position_in_parent_children->first))
      << "NewFork at nonexistent time " << time;

  // May be at |timeline_end()|.
  auto const timeline_it = timeline_find(time);

  // We cannot know the iterator into |this->children_| that the child object
  // will hold until after we have inserted it in |this->children_|.
  auto const child_it = children_.emplace(
      std::piecewise_construct,
      std::forward_as_tuple(time),
      std::forward_as_tuple(this /*parent*/,
                            children_.end(), /*position_in_parent_children*/
                            timeline_it /*position_in_parent_timeline*/));

  // Now set the iterator into |this->children_| in the child object.
  Forkable& child_forkable = child_it->second;
  child_forkable.position_in_parent_children = child_it;

  // Copy the tail of the trajectory in the child object.
  if (timeline_it != timeline_end()) {
    child_forkable.timeline_insert(++timeline_it, timeline_end());
  }
  return &child_forkable;
}

template<typename Tr4jectory>
void Forkable<Tr4jectory>::DeleteFork(not_null<Forkable**> const forkable) {
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

template<typename Tr4jectory>
bool Forkable<Tr4jectory>::is_root() const {
  return parent_ == nullptr;
}

template<typename Tr4jectory>
not_null<Forkable<Tr4jectory> const*> Forkable<Tr4jectory>::root() const {
  Forkable const* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Tr4jectory>
not_null<Forkable<Tr4jectory>*> Forkable<Tr4jectory>::root() {
  Forkable* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Tr4jectory>
Instant const* Forkable<Tr4jectory>::ForkTime() const {
  if (is_root()) {
    return nullptr;
  } else {
    return position_in_parent_children_->first;
  }
}

template<typename Tr4jectory>
Forkable<Tr4jectory>::Iterator Forkable<Tr4jectory>::Iterator::New(
    not_null<Forkable*> const forkable, Instant const & time) {
  Iterator iterator;

  // Go up the ancestry chain until we find a timeline that covers |time| (that
  // is, |time| is after the first time of the timeline).  Set |current_| to
  // the location of |time|, which may be |end()|.  The ancestry has |forkable|
  // at the back, and the object containing |current_| at the front.
  not_null<Forkable const*> ancestor = forkable;
  do {
    iterator.ancestry_.push_front(ancestor);
    if (!ancestor->timeline_is_empty() && ancestor->timeline_front() <= time) {
      iterator.current_ = ancestor->timeline_find(time);  // May be at end.
      break;
    }
    iterator.current_ = ancestor->timeline_end();
    ancestor = ancestor->parent_;
  } while (ancestor != nullptr);

  return iterator;
}

template<typename Tr4jectory>
Forkable<Tr4jectory>::Iterator Forkable<Tr4jectory>::Iterator::New(
    not_null<Forkable*> const forkable,
    not_null<Forkable*> const ancestor,
    typename Tr4jectory::TimelineConstIterator const
        position_in_ancestor_timeline) {
  Iterator iterator;

  // Go up the ancestry chain until we find |ancestor| and set |current_| to
  // |position_in_ancestor_timeline|.  The ancestry has |forkable|
  // at the back, and the object containing |current_| at the front.  If
  // |ancestor| is not found in the ancestry, stop at the root.
  not_null<Forkable const*> ancest0r = forkable;
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


template<typename Tr4jectory>
Forkable<Tr4jectory>::Iterator& Forkable<Tr4jectory>::Iterator::operator++() {
  CHECK(!ancestry_.empty());
  CHECK(current_ != ancestry_.front().end());

  // Check if there is a next child in the ancestry.
  auto ancestry_it = ++ancestry_.begin();
  if (ancestry_it != ancestry_.end()) {
    // There is a next child.  See if we reached its fork time.
    not_null<Forkable const*> child = *ancestry_it;
    Instant const& current_time = current_->first;
    Instant const& child_fork_time = child->position_in_ancestor_children;
    if (current_time == child_fork_time) {
      // Start iterating over the next child timeline.  Drop the leading
      // ancestor.
      current_ = child->timeline_first();
      ancestry_.pop_front();
      return *this;
    }
  }
  // Business as usual, keep moving along the same timeline.
  ++current_;

  return *this;
}

template<typename Tr4jectory>
typename Tr4jectory::TimelineConstIterator
Forkable<Tr4jectory>::current() const {
  return current_;
}

template<typename Tr4jectory>
Forkable<Tr4jectory>::Forkable(
    not_null<Forkable*> const parent,
    typename Children::const_iterator position_in_parent_children,
    typename Tr4jectory::TimelineConstIterator position_in_parent_timeline)
    : parent_(parent),
      position_in_parent_children_(position_in_parent_children),
      position_in_parent_timeline_(position_in_parent_timeline) {}

}  // namespace physics
}  // namespace principia
