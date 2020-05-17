
#pragma once

#include <deque>
#include <optional>
#include <vector>

#include "physics/forkable.hpp"

namespace principia {
namespace physics {
namespace internal_forkable {

template<typename Tr4jectory, typename It3rator, typename Traits>
not_null<Tr4jectory const*>
ForkableIterator<Tr4jectory, It3rator, Traits>::trajectory() const {
  CHECK(!ancestry_.empty());
  return ancestry_.back();
}

template<typename Tr4jectory, typename It3rator, typename Traits>
bool ForkableIterator<Tr4jectory, It3rator, Traits>::operator==(
    It3rator const& right) const {
  DCHECK_EQ(trajectory(), right.trajectory());
  // The comparison of iterators is faster than the comparison of deques, so if
  // this function returns false (which it does repeatedly in loops), it might
  // as well do so quickly.  There is a complication, however, because the two
  // iterators may not point to the same container, and we believe that
  // comparing them would be undefined behaviour; hence the size comparison,
  // which ensures that the two iterators are in the same fork and therefore can
  // legitimately be compared.
  return ancestry_.size() == right.ancestry_.size() &&
         current_ == right.current_ &&
         ancestry_ == right.ancestry_;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
bool ForkableIterator<Tr4jectory, It3rator, Traits>::operator!=(
    It3rator const& right) const {
  return !(*this == right);
}

template<typename Tr4jectory, typename It3rator, typename Traits>
It3rator& ForkableIterator<Tr4jectory, It3rator, Traits>::operator++() {
  CHECK(!ancestry_.empty());

  auto ancestry_it = ancestry_.begin();
  not_null<Tr4jectory const*> ancestor = *ancestry_it;
  auto current = ancestor->MakeEphemeral(current_);
  CHECK(current != ancestor->timeline_ephemeral_end());

  // Check if there is a next child in the ancestry.
  if (++ancestry_it != ancestry_.end()) {
    // There is a next child.  See if we reached its fork time.
    Instant const& current_time = Traits::time(current);
    not_null<Tr4jectory const*> child = *ancestry_it;
    Instant child_fork_time = (*child->position_in_parent_children_)->first;
    if (current_time == child_fork_time) {
      // We have reached the fork time of the next child.  There may be several
      // forks at that time so we must skip them until we find a fork that is at
      // a different time or the end of the children.
      do {
        current = child->timeline_ephemeral_begin();  // May be at end.
        ancestor = child;
        ancestry_.pop_front();
        if (++ancestry_it == ancestry_.end()) {
          break;
        }
        child = *ancestry_it;
        child_fork_time = (*child->position_in_parent_children_)->first;
      } while (current_time == child_fork_time);

      current_ = ancestor->MakeDurable(current);
      CheckNormalizedIfEnd();
      return *that();
    }
  }
  // Business as usual, keep moving along the same timeline.
  ++current;

  current_ = ancestor->MakeDurable(current);
  CheckNormalizedIfEnd();
  return *that();
}

template<typename Tr4jectory, typename It3rator, typename Traits>
It3rator& ForkableIterator<Tr4jectory, It3rator, Traits>::operator--() {
  CHECK(!ancestry_.empty());

  not_null<Tr4jectory const*> ancestor = ancestry_.front();
  auto current = ancestor->MakeEphemeral(current_);
  if (current == ancestor->timeline_ephemeral_begin()) {
    CHECK_NOTNULL(ancestor->parent_);
    // At the beginning of the first timeline.  Push the parent in front of the
    // ancestry and set |current_| to the fork point.  If the timeline is empty,
    // keep going until we find a non-empty one or the root.
    do {
      current = ancestor->parent_->MakeEphemeral(
          *ancestor->position_in_parent_timeline_);
      ancestor = ancestor->parent_;
      ancestry_.push_front(ancestor);
    } while (current == ancestor->timeline_ephemeral_end() &&
             ancestor->parent_ != nullptr);

    current_ = ancestor->MakeDurable(current);
    return *that();
  }
  --current;

  current_ = ancestor->MakeDurable(current);
  return *that();
}

template<typename Tr4jectory, typename It3rator, typename Traits>
typename ForkableIterator<Tr4jectory, It3rator, Traits>::
TimelineDurableConstIterator const&
ForkableIterator<Tr4jectory, It3rator, Traits>::current() const {
  return current_;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
void ForkableIterator<Tr4jectory, It3rator, Traits>::NormalizeIfEnd() {
  CHECK(!ancestry_.empty());
  if (current_ == ancestry_.front()->timeline_end() &&
      ancestry_.size() > 1) {
    ancestry_.erase(ancestry_.begin(), --ancestry_.end());
    current_ = ancestry_.front()->timeline_end();
  }
}

template<typename Tr4jectory, typename It3rator, typename Traits>
void ForkableIterator<Tr4jectory, It3rator, Traits>::CheckNormalizedIfEnd() {
  // Checking if the trajectory is a root is faster than obtaining the end of
  // the front of the deque, so it should be done first.
  CHECK(ancestry_.size() == 1 ||
        current_ != ancestry_.front()->timeline_end());
}

template<typename Tr4jectory, typename It3rator, typename Traits>
void Forkable<Tr4jectory, It3rator, Traits>::
DeleteFork(Tr4jectory*& trajectory) {
  CHECK_NOTNULL(trajectory);
  auto const fork_it = trajectory->Fork();
  // Find the position of |*trajectory| among our children and remove it.
  auto const range = children_.equal_range(Traits::time(fork_it.current_));
  for (auto it = range.first; it != range.second; ++it) {
    if (it->second.get() == trajectory) {
      children_.erase(it);
      trajectory = nullptr;
      return;
    }
  }
  LOG(FATAL) << "argument is not a child of this trajectory";
}

template<typename Tr4jectory, typename It3rator, typename Traits>
bool Forkable<Tr4jectory, It3rator, Traits>::is_root() const {
  return parent_ == nullptr;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
not_null<Tr4jectory const*>
Forkable<Tr4jectory, It3rator, Traits>::root() const {
  not_null<Tr4jectory const*> ancestor = that();
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
not_null<Tr4jectory*> Forkable<Tr4jectory, It3rator, Traits>::root() {
  not_null<Tr4jectory*> ancestor = that();
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
not_null<Tr4jectory const*> Forkable<Tr4jectory, It3rator, Traits>::
parent() const {
  return parent_;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
not_null<Tr4jectory*> Forkable<Tr4jectory, It3rator, Traits>::parent() {
  return parent_;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
It3rator Forkable<Tr4jectory, It3rator, Traits>::begin() const {
  not_null<Tr4jectory const*> ancestor = root();
  return Wrap(ancestor, ancestor->timeline_begin());
}

template<typename Tr4jectory, typename It3rator, typename Traits>
It3rator Forkable<Tr4jectory, It3rator, Traits>::end() const {
  not_null<Tr4jectory const*> const ancestor = that();
  It3rator iterator;
  iterator.ancestry_.push_front(ancestor);
  iterator.current_ = ancestor->timeline_end();
  iterator.CheckNormalizedIfEnd();
  return iterator;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
typename It3rator::reference
Forkable<Tr4jectory, It3rator, Traits>::front() const {
  // TODO(phl): This can be implemented more efficiently.
  return *begin();
}

template<typename Tr4jectory, typename It3rator, typename Traits>
typename It3rator::reference Forkable<Tr4jectory, It3rator, Traits>::
back() const {
  // TODO(phl): This can be implemented more efficiently.
  return *--end();
}

template<typename Tr4jectory, typename It3rator, typename Traits>
It3rator Forkable<Tr4jectory, It3rator, Traits>::
Find(Instant const& time) const {
  It3rator iterator;

  // Go up the ancestry chain until we find a timeline that covers |time| (that
  // is, |time| is after the first time of the timeline).  Set |current_| to
  // the location of |time|, which may be |end()|.  The ancestry has |forkable|
  // at the back, and the object containing |current_| at the front.
  Tr4jectory const* ancestor = that();
  TimelineEphemeralConstIterator current;
  do {
    iterator.ancestry_.push_front(ancestor);
    if (!ancestor->timeline_empty() &&
        Traits::time(ancestor->timeline_ephemeral_begin()) <= time) {
      current =
          ancestor->timeline_ephemeral_find(time);  // May be at end.
      break;
    }
    current = ancestor->timeline_ephemeral_end();
    ancestor = ancestor->parent_;
  } while (ancestor != nullptr);

  iterator.current_ = iterator.ancestry_.front()->MakeDurable(current);
  iterator.NormalizeIfEnd();
  return iterator;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
It3rator Forkable<Tr4jectory, It3rator, Traits>::
LowerBound(Instant const& time) const {
  It3rator iterator;
  Tr4jectory const* ancestor = that();

  // This queue is parallel to |iterator.ancestry_|, and each entry is an
  // iterator in the corresponding ancestry timeline.  Note that we use a
  // |nullopt| sentinel for the innermost timeline.
  std::deque<std::optional<TimelineDurableConstIterator>> fork_points;
  fork_points.push_front(std::nullopt);

  // Go up the ancestry chain until we find a (nonempty) timeline that covers
  // |time| (that is, |time| is on or after the first time of the timeline).
  TimelineEphemeralConstIterator current;
  do {
    iterator.ancestry_.push_front(ancestor);
    if (!ancestor->timeline_empty() &&
        Traits::time(ancestor->timeline_ephemeral_begin()) <= time) {
      // We have found a timeline that covers |time|.  Find where |time| falls
      // in that timeline (that may be after the end).

      current = ancestor->timeline_ephemeral_lower_bound(time);

      // Check if the returned iterator is directly usable.
      auto const& fork_point = fork_points.front();
      if (current == ancestor->timeline_ephemeral_end() ||
          (fork_point &&
           *fork_point != ancestor->timeline_end() &&
           Traits::time(*fork_point) < Traits::time(current))) {
        // |time| is after the end of this timeline or after the |fork_point|
        // (if any).  We may have to return an |End| iterator, so let's prepare
        // |current| for that case.
        current = ancestor->timeline_ephemeral_end();

        // Check if we have a more nested fork with a point before |time|.  Go
        // down the ancestry looking for a timeline that is nonempty and not
        // forked at the same point as its parent.
        auto ancestry_it = iterator.ancestry_.begin();
        auto fork_points_it = fork_points.begin();
        for (;;) {
          ++ancestry_it;
          ++fork_points_it;
          if (ancestry_it == iterator.ancestry_.end()) {
            // We didn't find an interesting fork in the ancestry, so we stop
            // here and |NormalizeIfEnd| will return a proper |End|.
            CHECK(fork_points_it == fork_points.end());
            break;
          }
          if (!(*ancestry_it)->timeline_empty() &&
                (!*fork_points_it ||
                 **fork_points_it != (*ancestry_it)->timeline_end())) {
            // We found an interesting timeline, i.e. one that is nonempty and
            // not forked at the fork point of its parent.  Cut the ancestry and
            // return the beginning of that timeline.
            iterator.ancestry_.erase(iterator.ancestry_.begin(), ancestry_it);
            current = (*ancestry_it)->timeline_ephemeral_begin();
            break;
          }
        }
      }
      break;
    }
    fork_points.push_front(ancestor->position_in_parent_timeline_);
    current = ancestor->timeline_ephemeral_begin();
    ancestor = ancestor->parent_;
  } while (ancestor != nullptr);

  iterator.current_ = iterator.ancestry_.front()->MakeDurable(current);
  iterator.NormalizeIfEnd();
  return iterator;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
It3rator Forkable<Tr4jectory, It3rator, Traits>::Fork() const {
  CHECK(!is_root());
  not_null<Tr4jectory const*> ancestor = that();
  TimelineDurableConstIterator position_in_ancestor_timeline;
  do {
    position_in_ancestor_timeline = *ancestor->position_in_parent_timeline_;
    ancestor = ancestor->parent_;
  } while (position_in_ancestor_timeline == ancestor->timeline_end() &&
           ancestor->parent_ != nullptr);
  return Wrap(ancestor, position_in_ancestor_timeline);
}

template<typename Tr4jectory, typename It3rator, typename Traits>
std::int64_t Forkable<Tr4jectory, It3rator, Traits>::Size() const {
  Tr4jectory const* ancestor = that();

  // Get the size directly for the leaf trajectory, this is more efficient if
  // there are no forks.
  std::int64_t size = ancestor->timeline_size();

  // Go up the ancestry chain adding the sizes.
  Tr4jectory const* parent = ancestor->parent_;
  while (parent != nullptr) {
    if (!parent->timeline_empty()) {
      size += std::distance(parent->timeline_ephemeral_begin(),
                            ancestor->MakeEphemeral(
                                *ancestor->position_in_parent_timeline_)) + 1;
    }
    ancestor = parent;
    parent = ancestor->parent_;
  }

  return size;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
bool Forkable<Tr4jectory, It3rator, Traits>::Empty() const {
  // If this object has an ancestor surely it is hooked off of a point in some
  // timeline, so this object is not empty.
  return timeline_empty() && parent_ == nullptr;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
not_null<Tr4jectory*> Forkable<Tr4jectory, It3rator, Traits>::NewFork(
    TimelineEphemeralConstIterator const& timeline_it) {
  // First create a child in the multimap.
  Instant time;
  if (timeline_it == timeline_ephemeral_end()) {
    CHECK(!is_root());
    time = (*position_in_parent_children_)->first;
  } else {
    time = Traits::time(timeline_it);
  }
  auto const child_it = children_.emplace(time, std::make_unique<Tr4jectory>());

  // Now set the members of the child object.
  std::unique_ptr<Tr4jectory> const& child_forkable = child_it->second;
  child_forkable->parent_ = that();
  child_forkable->position_in_parent_children_ = child_it;
  child_forkable->position_in_parent_timeline_ = MakeDurable(timeline_it);

  return child_forkable.get();
}

template<typename Tr4jectory, typename It3rator, typename Traits>
void Forkable<Tr4jectory, It3rator, Traits>::AttachForkToCopiedBegin(
    not_null<std::unique_ptr<Tr4jectory>> fork) {
  CHECK(fork->is_root());
  CHECK(!fork->timeline_empty());
  auto const fork_timeline_begin = fork->timeline_begin();
  auto const fork_timeline_end = fork->timeline_end();

  // The children of |fork| whose |position_in_parent_timeline_| was at
  // |begin()| are referencing a point that will soon be removed from the
  // timeline.  They must now point at |end()| to indicate that their fork time
  // is not in |fork|'s timeline.
  for (auto const& [_, child] : fork->children_) {
    if (child->position_in_parent_timeline_ == fork_timeline_begin) {
      child->position_in_parent_timeline_ = fork_timeline_end;
    }
  }

  // Insert |fork| in the |children_| of this object.
  auto const child_it = children_.emplace_hint(
      children_.end(),
      Traits::time(fork_timeline_begin),
      std::move(fork));

  // Set the pointer into this object.  Note that |fork| is no longer usable.
  child_it->second->parent_ = that();
  child_it->second->position_in_parent_children_ = child_it;
  auto position_in_parent_timeline = timeline_ephemeral_end();
  if (!timeline_empty()) {
    --position_in_parent_timeline;
  }
  child_it->second->position_in_parent_timeline_ =
      MakeDurable(position_in_parent_timeline);
}

template<typename Tr4jectory, typename It3rator, typename Traits>
not_null<std::unique_ptr<Tr4jectory>>
Forkable<Tr4jectory, It3rator, Traits>::DetachForkWithCopiedBegin() {
  CHECK(!is_root());

  // The children whose |position_in_parent_timeline_| was at |end()| are those
  // whose fork time was not in this object's timeline.  The caller must have
  // ensured that now it is, so point them to the beginning of this timeline.
  for (auto const& [_, child] : children_) {
    if (child->position_in_parent_timeline_ == timeline_end()) {
      child->position_in_parent_timeline_ = timeline_begin();
    }
  }

  // Remove this trajectory from the children of its parent.
  auto owned_this = std::move((*position_in_parent_children_)->second);
  parent_->children_.erase(*position_in_parent_children_);

  // Clear all the pointers to the parent.
  parent_ = nullptr;
  position_in_parent_children_ = std::nullopt;
  position_in_parent_timeline_ = std::nullopt;

  return std::move(owned_this);
}

template<typename Tr4jectory, typename It3rator, typename Traits>
void Forkable<Tr4jectory, It3rator, Traits>::
DeleteAllForksAfter(Instant const& time) {
  // Get an iterator denoting the first entry with time > |time|.  Remove that
  // entry and all the entries that follow it.  This preserves any entry with
  // time == |time|.
  CHECK(is_root() || time >= Traits::time(Fork().current_))
      << "DeleteAllForksAfter before the fork time " << time << " vs. "
      << Traits::time(Fork().current_);
  auto const it = children_.upper_bound(time);
  children_.erase(it, children_.end());
}

template<typename Tr4jectory, typename It3rator, typename Traits>
void Forkable<Tr4jectory, It3rator, Traits>::
CheckNoForksBefore(Instant const& time) {
  CHECK(is_root()) << "CheckNoForksBefore on a nonroot trajectory";
  // Get an iterator denoting the first entry with time >= |time|.  Check that
  // there are no forks before it.  A fork with time == |time| is fine.
  auto const it = children_.lower_bound(time);
  CHECK(children_.begin() == it) << "CheckNoForksBefore found "
                                 << std::distance(children_.begin(), it)
                                 << " forks before " << time;
}

template<typename Tr4jectory, typename It3rator, typename Traits>
void Forkable<Tr4jectory, It3rator, Traits>::WriteSubTreeToMessage(
    not_null<serialization::DiscreteTrajectory*> const message,
    std::vector<Tr4jectory*>& forks) const {
  std::optional<Instant> last_instant;
  serialization::DiscreteTrajectory::Litter* litter = nullptr;
  for (auto const& [fork_time, child] : children_) {
    // Determine if this |child| needs to be serialized.  If so, record its
    // position in |fork_positions| and null out its pointer in |forks|.
    // Apologies for the O(N) search.
    auto const it = std::find(forks.begin(), forks.end(), child.get());
    if (it == forks.end()) {
      continue;
    } else {
      message->add_fork_position(it - forks.begin());
      *it = nullptr;
    }

    if (!last_instant || fork_time != last_instant) {
      last_instant = fork_time;
      litter = message->add_children();
      fork_time.WriteToMessage(litter->mutable_fork_time());
    }
    child->WriteSubTreeToMessage(litter->add_trajectories(), forks);
  }
}

template<typename Tr4jectory, typename It3rator, typename Traits>
void Forkable<Tr4jectory, It3rator, Traits>::FillSubTreeFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<Tr4jectory**> const& forks) {
  // There were no fork positions prior to Буняковский.
  bool const has_fork_position = message.fork_position_size() > 0;
  std::int32_t index = 0;
  for (serialization::DiscreteTrajectory::Litter const& litter :
           message.children()) {
    Instant const fork_time = Instant::ReadFromMessage(litter.fork_time());
    for (serialization::DiscreteTrajectory const& child :
             litter.trajectories()) {
      not_null<Tr4jectory*> fork = NewFork(timeline_ephemeral_find(fork_time));
      fork->FillSubTreeFromMessage(child, forks);
      if (has_fork_position) {
        std::int32_t const fork_position = message.fork_position(index);
        *forks[fork_position] = fork;
      }
      ++index;
    }
  }
}

template<typename Tr4jectory, typename It3rator, typename Traits>
It3rator Forkable<Tr4jectory, It3rator, Traits>::Wrap(
    not_null<const Tr4jectory*> const ancestor,
    TimelineDurableConstIterator const position_in_ancestor_timeline) const {
  It3rator iterator;

  // Go up the ancestry chain until we find |ancestor| and set |current_| to
  // |position_in_ancestor_timeline|.  The ancestry has |that()| at the back,
  // and the object containing |current_| at the front.
  not_null<Tr4jectory const*> ancest0r = that();
  do {
    iterator.ancestry_.push_front(ancest0r);
    if (ancestor == ancest0r) {
      iterator.current_ = position_in_ancestor_timeline;  // May be at end.
      iterator.CheckNormalizedIfEnd();
      return iterator;
    }
    ancest0r = ancest0r->parent_;
  } while (ancest0r != nullptr);

  LOG(FATAL) << "The ancestor parameter is not an ancestor of this trajectory";
  base::noreturn();
}

}  // namespace internal_forkable
}  // namespace physics
}  // namespace principia
