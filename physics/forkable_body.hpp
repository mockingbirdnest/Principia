#include "forkable.hpp"
#pragma once

#include "physics/forkable.hpp"

namespace principia {
namespace physics {

template<typename Tr4jectory>
not_null<Forkable<Tr4jectory>*> Forkable<Tr4jectory>::NewFork(
    Instant const & time) {
  CHECK(timeline_.find(time) != timeline_.end() ||
        (!is_root() && time == ForkTime()))
      << "NewFork at nonexistent time " << time;

  // May be at |end()|.
  auto fork_it = timeline_.find(time);

  // We cannot know the iterator into children_ until after we have done the
  // insertion in children_.
  Fork const fork = {children_.end(), fork_it};
  auto const child_it = children_.emplace(
      std::piecewise_construct,
      std::forward_as_tuple(time),
      std::forward_as_tuple(body_, this /*parent*/, fork));
  if (fork_it != timeline_.end()) {
    child_it->second.timeline_.insert(++fork_it, timeline_.end());
  }
  child_it->second.fork_->children = child_it;
  return &child_it->second;
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
Instant const* Forkable<Tr4jectory>::fork_time() const {
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

}  // namespace physics
}  // namespace principia
