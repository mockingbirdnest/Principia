#include "forkable.hpp"
#pragma once

#include "physics/forkable.hpp"

namespace principia {
namespace physics {

template<typename Tr4jectory>
not_null<Forkable<Tr4jectory>*> Forkable<Tr4jectory>::NewFork(
    Instant const & time) {
  CHECK(ContainsTime(time)) << "NewFork at nonexistent time " << time;

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
void Forkable<Tr4jectory>::DeleteFork(not_null<Forkable**> const fork) {
  CHECK_NOTNULL(*fork);
  Instant const* const fork_time = (*fork)->fork_time();
  CHECK_NOTNULL(fork_time);
  // Find the position of |*fork| among our children and remove it.
  auto const range = children_.equal_range(*fork_time);
  for (auto it = range.first; it != range.second; ++it) {
    if (&it->second == *fork) {
      children_.erase(it);
      *fork = nullptr;
      return;
    }
  }
  LOG(FATAL) << "fork is not a child of this trajectory";
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
bool Forkable<Tr4jectory>::ContainsTime(Instant const & time) const {
  if (timeline_find(time) != timeline_end()) {
    return true
  } else if (is_root()) {
    return false;
  } else {
    return time == *ForkTime();
  }
}

template<typename Tr4jectory>
Instant const* Forkable<Tr4jectory>::ForkTime() const {
  not_null<Forkable const*> ancestor = this;
  while (ancestor->parent_ != nullptr) {
    if (!ancestor->parent_->timeline_is_empty()) {
      return &(ancestor->position_in_parent_timeline_->first);
    }
    ancestor = ancestor->parent_;
  }
  return nullptr;
}

template<typename Tr4jectory>
Iterator Forkable<Tr4jectory>::Iterator::New(
    not_null<Forkable*> const forkable, Instant const & time) {
  not_null<Forkable const*> ancestor = forkable;
  while (ancestor->parent_ != nullptr &&
         (ancestor->fork_->timeline == ancestor->parent_->timeline_.end() ||
          time <= ancestor->fork_->timeline->first)) {
    ancestry_.push_front(ancestor);
    forks_.push_front(*ancestor->fork_);
    ancestor = ancestor->parent_;
  }
  ancestry_.push_front(ancestor);
  current_ = ancestor->timeline_.lower_bound(time);
  CHECK(!current_is_misplaced());
}

template<typename Tr4jectory>
Iterator Forkable<Tr4jectory>::Iterator::New(
    not_null<Forkable*> const forkable,
    not_null<Forkable*> const ancestor,
    typename Tr4jectory::TimelineConstIterator const
        position_in_ancestor_timeline) {
  return Iterator();
}

template<typename Tr4jectory>
Forkable<Tr4jectory>::Forkable(
    not_null<Forkable*> const parent,
    typename Children::const_iterator position_in_parent_children,
    typename Tr4jectory::TimelineConstIterator position_in_parent_timeline) {
}

}  // namespace physics
}  // namespace principia
