#pragma once

#include "physics/forkable.hpp"

namespace principia {
namespace physics {

template<typename Tr4jectory>
not_null<Tr4jectory*> Forkable<Tr4jectory>::NewFork(Instant const & time) {
  Instant const* const fork_time = ForkTime();
  CHECK(timeline_find(time) != timeline_end() ||
        (fork_time != nullptr && time == *fork_time))
      << "NewFork at nonexistent time " << time;

  // May be at |timeline_end()| if |time| is the fork time of this object.
  auto timeline_it = timeline_find(time);

  // First create a child in the multimap.
  auto const child_it = children_.emplace(std::piecewise_construct,
                                          std::forward_as_tuple(time),
                                          std::forward_as_tuple());

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

template<typename Tr4jectory>
void Forkable<Tr4jectory>::DeleteFork(not_null<Tr4jectory**> const trajectory) {
  CHECK_NOTNULL(*trajectory);
  Instant const* const fork_time = (*trajectory)->ForkTime();
  CHECK_NOTNULL(fork_time);
  // Find the position of |*forkable| among our children and remove it.
  auto const range = children_.equal_range(*fork_time);
  for (auto it = range.first; it != range.second; ++it) {
    if (&it->second == *trajectory) {
      children_.erase(it);
      *trajectory = nullptr;
      return;
    }
  }
  LOG(FATAL) << "argument is not a child of this trajectory";
}

template<typename Tr4jectory>
void Forkable<Tr4jectory>::ForgetAfter(Instant const& time) {
  // Each of these blocks gets an iterator denoting the first entry with
  // time > |time|.  It then removes that entry and all the entries that follow
  // it.  This preserve any entry with time == |time|.
  {
    auto const it = timeline_upper_bound(time);
    CHECK(is_root() || time >= *ForkTime())
        << "ForgetAfter before the fork time";
    timeline_erase(it, timeline_end());
  }
  {
    auto const it = children_.upper_bound(time);
    children_.erase(it, children_.end());
  }
}

template<typename Tr4jectory>
void Forkable<Tr4jectory>::ForgetBefore(Instant const& time) {
  // Check that this is a root.
  CHECK(is_root()) << "ForgetBefore on a nonroot trajectory";
  // Each of these blocks gets an iterator denoting the first entry with
  // time > |time|.  It then removes that all the entries that precede it.  This
  // removes any entry with time == |time|.
  {
    auto it = timeline_upper_bound(time);
    timeline_erase(timeline_begin(), it);
  }
  {
    auto it = children_.upper_bound(time);
    children_.erase(children_.begin(), it);
  }
}
template<typename Tr4jectory>
bool Forkable<Tr4jectory>::is_root() const {
  return parent_ == nullptr;
}

template<typename Tr4jectory>
not_null<Tr4jectory const*>
Forkable<Tr4jectory>::root() const {
  not_null<Tr4jectory const*> ancestor = that();
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Tr4jectory>
not_null<Tr4jectory*>
Forkable<Tr4jectory>::root() {
  not_null<Tr4jectory*> ancestor = that();
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
    return &position_in_parent_children_->first;
  }
}

template<typename Tr4jectory>
bool Forkable<Tr4jectory>::Iterator::operator==(Iterator const& right) const {
  return ancestry_ == right.ancestry_ && current_ == right.current_;
}

template<typename Tr4jectory>
bool Forkable<Tr4jectory>::Iterator::operator!=(Iterator const& right) const {
  return !(*this == right);
}

template<typename Tr4jectory>
typename Forkable<Tr4jectory>::Iterator&
Forkable<Tr4jectory>::Iterator::operator++() {
  CHECK(!ancestry_.empty());
  CHECK(current_ != ancestry_.front()->timeline_end());

  // Check if there is a next child in the ancestry.
  auto ancestry_it = ancestry_.begin();
  if (++ancestry_it != ancestry_.end()) {
    // There is a next child.  See if we reached its fork time.
    Instant const& current_time = ForkableTraits<Tr4jectory>::time(current_);
    not_null<Tr4jectory const*> child = *ancestry_it;
    Instant child_fork_time = child->position_in_parent_children_->first;
    if (current_time == child_fork_time) {
      // We have reached the fork time of the next child.  There may be several
      // forks at that time so we must skip them until we find a fork that is at
      // a different time or the end of the children.
      do {
        current_ = child->timeline_begin();  // May be at end.
        ancestry_.pop_front();
        if (++ancestry_it == ancestry_.end()) {
          break;
        }
        child = *ancestry_it;
        child_fork_time = child->position_in_parent_children_->first;
      } while (current_time == child_fork_time);

      CheckNormalizedIfEnd();
      return *this;
    }
  }
  // Business as usual, keep moving along the same timeline.
  ++current_;

  CheckNormalizedIfEnd();
  return *this;
}

template<typename Tr4jectory>
typename Forkable<Tr4jectory>::Iterator&
Forkable<Tr4jectory>::Iterator::operator--() {
  CHECK(!ancestry_.empty());

  not_null<Tr4jectory const*> ancestor = ancestry_.front();
  if (current_ == ancestor->timeline_begin()) {
    CHECK_NOTNULL(ancestor->parent_);
    // At the beginning of the first timeline.  Push the parent in front of the
    // ancestry and set |current_| to the fork point.  If the timeline is empty,
    // keep going until we find a non-empty one or the root.
    do {
      current_ = ancestor->position_in_parent_timeline_;
      ancestor = ancestor->parent_;
      ancestry_.push_front(ancestor);
    } while (current_ == ancestor->timeline_end() &&
             ancestor->parent_ != nullptr);
    return *this;
  }

  --current_;
  return *this;
}

template<typename Tr4jectory>
typename Forkable<Tr4jectory>::TimelineConstIterator
Forkable<Tr4jectory>::Iterator::current() const {
  return current_;
}

template<typename Tr4jectory>
not_null<Tr4jectory const*> Forkable<Tr4jectory>::Iterator::trajectory() const {
  CHECK(!ancestry_.empty());
  return ancestry_.back();
}

template<typename Tr4jectory>
void Forkable<Tr4jectory>::Iterator::NormalizeIfEnd() {
  CHECK(!ancestry_.empty());
  if (current_ == ancestry_.front()->timeline_end() &&
      ancestry_.size() > 1) {
    ancestry_.erase(ancestry_.begin(), --ancestry_.end());
    current_ = ancestry_.front()->timeline_end();
  }
}

template<typename Tr4jectory>
void Forkable<Tr4jectory>::Iterator::CheckNormalizedIfEnd() {
  CHECK(current_ != ancestry_.front()->timeline_end() ||
        ancestry_.size() == 1);
}

template<typename Tr4jectory>
typename Forkable<Tr4jectory>::Iterator
Forkable<Tr4jectory>::Begin() const {
  not_null<Tr4jectory const*> ancestor = root();
  return Wrap(ancestor, ancestor->timeline_begin());
}

template<typename Tr4jectory>
typename Forkable<Tr4jectory>::Iterator
Forkable<Tr4jectory>::End() const {
  not_null<Tr4jectory const*> ancestor = that();
  return Wrap(ancestor, ancestor->timeline_end());
}

template<typename Tr4jectory>
typename Forkable<Tr4jectory>::Iterator
Forkable<Tr4jectory>::Find(Instant const& time) const {
  Iterator iterator;

  // Go up the ancestry chain until we find a timeline that covers |time| (that
  // is, |time| is after the first time of the timeline).  Set |current_| to
  // the location of |time|, which may be |end()|.  The ancestry has |forkable|
  // at the back, and the object containing |current_| at the front.
  Tr4jectory const* ancestor = that();
  do {
    iterator.ancestry_.push_front(ancestor);
    if (!ancestor->timeline_empty() &&
        ForkableTraits<Tr4jectory>::time(ancestor->timeline_begin()) <= time) {
      iterator.current_ = ancestor->timeline_find(time);  // May be at end.
      break;
    }
    iterator.current_ = ancestor->timeline_end();
    ancestor = ancestor->parent_;
  } while (ancestor != nullptr);

  iterator.NormalizeIfEnd();
  return iterator;
}

template<typename Tr4jectory>
typename Forkable<Tr4jectory>::Iterator
Forkable<Tr4jectory>::LowerBound(Instant const& time) const {
  Iterator iterator;

  // Go up the ancestry chain until we find a timeline that covers |time| (that
  // is, |time| is after the first time of the timeline).  Set |current_| to
  // the location of |time|, which may be |end()|.  The ancestry has |forkable|
  // at the back, and the object containing |current_| at the front.
  Tr4jectory const* ancestor = that();
  do {
    iterator.ancestry_.push_front(ancestor);
    if (!ancestor->timeline_empty() &&
        ForkableTraits<Tr4jectory>::time(ancestor->timeline_begin()) <= time) {
      iterator.current_ = ancestor->timeline_lower_bound(time);  // May be at end.
      break;
    }
    iterator.current_ = ancestor->timeline_begin();
    ancestor = ancestor->parent_;
  } while (ancestor != nullptr);

  iterator.NormalizeIfEnd();
  return iterator;
}

template<typename Tr4jectory>
typename Forkable<Tr4jectory>::Iterator
Forkable<Tr4jectory>::Wrap(
    not_null<const Tr4jectory*> const ancestor,
    TimelineConstIterator const position_in_ancestor_timeline) const {
  Iterator iterator;

  // Go up the ancestry chain until we find |ancestor| and set |current_| to
  // |position_in_ancestor_timeline|.  The ancestry has |forkable|
  // at the back, and the object containing |current_| at the front.
  not_null<Tr4jectory const*> ancest0r = that();
  do {
    iterator.ancestry_.push_front(ancest0r);
    if (ancestor == ancest0r) {
      iterator.current_ = position_in_ancestor_timeline;  // May be at end.
      iterator.CheckNormalizedIfEnd();
      return iterator;
    }
    iterator.current_ = ancest0r->timeline_end();
    ancest0r = ancest0r->parent_;
  } while (ancest0r != nullptr);

  LOG(FATAL) << "The ancestor parameter is not an ancestor of this trajectory";
  base::noreturn();
}

template<typename Tr4jectory>
void Forkable<Tr4jectory>::WritePointerToMessage(
    not_null<serialization::Trajectory::Pointer*> const message) const {
  not_null<Tr4jectory const*> ancestor = that();
  while (ancestor->parent_ != nullptr) {
    auto const position_in_parent_children = position_in_parent_children_;
    auto const position_in_parent_timeline = position_in_parent_timeline_;
    ancestor = ancestor->parent_;
    int const children_distance =
        std::distance(ancestor->children_.begin(), position_in_parent_children);
    int const timeline_distance =
        std::distance(ancestor->timeline_begin(), position_in_parent_timeline);
    auto* const fork_message = message->add_fork();
    fork_message->set_children_distance(children_distance);
    fork_message->set_timeline_distance(timeline_distance);
  }
}

template<typename Tr4jectory>
not_null<Tr4jectory*> Forkable<Tr4jectory>::ReadPointerFromMessage(
    serialization::Trajectory::Pointer const& message,
    not_null<Tr4jectory*> const trajectory) {
  CHECK(trajectory->is_root());
  not_null<Tr4jectory*> descendant = trajectory;
  for (auto const& fork_message : message.fork()) {
    int const children_distance = fork_message.children_distance();
    int const timeline_distance = fork_message.timeline_distance();
    auto children_it = descendant->children_.begin();
    auto timeline_it = descendant->timeline_begin();
    std::advance(children_it, children_distance);
    std::advance(timeline_it, timeline_distance);
    descendant = &children_it->second;
  }
  return descendant;
}

template<typename Tr4jectory>
void Forkable<Tr4jectory>::WriteSubTreeToMessage(
    not_null<serialization::Trajectory*> const message) const {
  Instant last_instant;  // optional.
  bool is_first = true;
  serialization::Trajectory::Litter* litter = nullptr;
  for (auto const& pair : children_) {
    Instant const& fork_time = pair.first;
    Tr4jectory const& child = pair.second;
    if (is_first || fork_time != last_instant) {
      is_first = false;
      last_instant = fork_time;
      litter = message->add_children();
      fork_time.WriteToMessage(litter->mutable_fork_time());
    }
    child.WriteSubTreeToMessage(litter->add_trajectories());
  }
}

}  // namespace physics
}  // namespace principia
