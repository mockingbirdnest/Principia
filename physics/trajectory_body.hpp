#pragma once

#include "trajectory.hpp"

#include <algorithm>
#include <list>
#include <map>

#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "physics/oblate_body.hpp"

namespace principia {

using base::make_not_null_unique;
using geometry::Instant;

namespace physics {

template<typename Frame>
Trajectory<Frame>::Trajectory(not_null<Body const*> const body)
    : body_(body),
      parent_(nullptr) {
  CHECK(body_->is_compatible_with<Frame>())
      << "Oblate body not in the same frame as the trajectory";
}

template<typename Frame>
typename Trajectory<Frame>::NativeIterator Trajectory<Frame>::first() const {
  NativeIterator it;
  it.InitializeFirst(this);
  return it;
}

template<typename Frame>
typename Trajectory<Frame>::NativeIterator Trajectory<Frame>::on_or_after(
    Instant const& time) const {
  NativeIterator it;
  it.InitializeOnOrAfter(time, this);
  return it;
}

template<typename Frame>
typename Trajectory<Frame>::NativeIterator Trajectory<Frame>::last() const {
  NativeIterator it;
  it.InitializeLast(this);
  return it;
}

template<typename Frame>
template<typename ToFrame>
typename Trajectory<Frame>::TEMPLATE TransformingIterator<ToFrame>
Trajectory<Frame>::first_with_transform(
    Transform<ToFrame> const& transform) const {
  TransformingIterator<ToFrame> it(transform);
  it.InitializeFirst(this);
  return it;
}

template<typename Frame>
template<typename ToFrame>
typename Trajectory<Frame>::TEMPLATE TransformingIterator<ToFrame>
Trajectory<Frame>::on_or_after_with_transform(
    Instant const& time,
    Transform<ToFrame> const& transform) const {
  TransformingIterator<ToFrame> it(transform);
  it.InitializeOnOrAfter(time, this);
  return it;
}

template<typename Frame>
template<typename ToFrame>
typename Trajectory<Frame>::TEMPLATE TransformingIterator<ToFrame>
Trajectory<Frame>::last_with_transform(
    Transform<ToFrame> const& transform) const {
  TransformingIterator<ToFrame> it(transform);
  it.InitializeLast(this);
  return it;
}

template<typename Frame>
std::map<Instant, Position<Frame>> Trajectory<Frame>::Positions() const {
  std::map<Instant, Position<Frame>> result;
  for (NativeIterator it = first(); !it.at_end(); ++it) {
    Instant const& time = it.time();
    result.emplace_hint(result.end(), time, it.degrees_of_freedom().position());
  }
  return result;
}

template<typename Frame>
std::map<Instant, Velocity<Frame>> Trajectory<Frame>::Velocities() const {
  std::map<Instant, Velocity<Frame>> result;
  for (NativeIterator it = first(); !it.at_end(); ++it) {
    Instant const& time = it.time();
    result.emplace_hint(result.end(), time, it.degrees_of_freedom().velocity());
  }
  return result;
}

template<typename Frame>
std::list<Instant> Trajectory<Frame>::Times() const {
  std::list<Instant> result;
  for (NativeIterator it = first(); !it.at_end(); ++it) {
    result.push_back(it.time());
  }
  return result;
}

template<typename Frame>
void Trajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  auto inserted = timeline_.emplace(time, degrees_of_freedom);
  CHECK(timeline_.end() == ++inserted.first) << "Append out of order";
  CHECK(inserted.second) << "Append at existing time " << time
                         << ", time range = [" << Times().front() << ", "
                         << Times().back() << "]";
}

template<typename Frame>
void Trajectory<Frame>::ForgetAfter(Instant const& time) {
  CHECK(InTimelineOrAtFork(time)) << "ForgetAfter a nonexistent time " << time;

  // Each of these blocks gets an iterator denoting the first entry with
  // time > |time|.  It then removes that entry and all the entries that follow
  // it.  This preserve any entry with time == |time|.
  {
    auto const it = timeline_.upper_bound(time);
    timeline_.erase(it, timeline_.end());
  }
  {
    auto const it = children_.upper_bound(time);
    children_.erase(it, children_.end());
  }
}

template<typename Frame>
void Trajectory<Frame>::ForgetBefore(Instant const& time) {
  // Check that this is a root.
  CHECK(is_root()) << "ForgetBefore on a nonroot trajectory";
  // Check that |time| is the time of one of our Timeline.
  CHECK(timeline_.find(time) != timeline_.end())
      << "ForgetBefore a nonexistent time";
  {
    auto it = timeline_.upper_bound(time);
    timeline_.erase(timeline_.begin(), it);
  }
  {
    auto it = children_.upper_bound(time);
    children_.erase(children_.begin(), it);
  }
}

template<typename Frame>
not_null<Trajectory<Frame>*> Trajectory<Frame>::NewFork(Instant const& time) {
  CHECK(InTimelineOrAtFork(time)) << "NewFork at nonexistent time " << time;

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

template<typename Frame>
void Trajectory<Frame>::DeleteFork(not_null<Trajectory**> const fork) {
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


template<typename Frame>
bool Trajectory<Frame>::is_root() const {
  return parent_ == nullptr;
}

template<typename Frame>
not_null<Trajectory<Frame> const*> Trajectory<Frame>::root() const {
  Trajectory const* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Frame>
not_null<Trajectory<Frame>*> Trajectory<Frame>::root() {
  Trajectory* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Frame>
Instant const* Trajectory<Frame>::fork_time() const {
  if (parent_ == nullptr) {
    return nullptr;
  } else {
    return &(fork_->timeline->first);
  }
}

template<typename Frame>
template<typename B>
std::enable_if_t<std::is_base_of<Body, B>::value, not_null<B const*>>
Trajectory<Frame>::body() const {
// Dynamic casting is expensive, as in 3x slower for the benchmarks.  Do that in
// debug mode to catch bugs, but not in optimized mode where we want all the
// performance we can get.
#ifdef _DEBUG
  return dynamic_cast<B const*>(static_cast<Body const*>(body_));
#else
  return static_cast<not_null<B const*>>(body_);
#endif
}

template<typename Frame>
void Trajectory<Frame>::set_intrinsic_acceleration(
    IntrinsicAcceleration const acceleration) {
  CHECK(body_->is_massless()) << "Trajectory is for a massive body";
  CHECK(intrinsic_acceleration_ == nullptr)
      << "Trajectory already has an intrinsic acceleration";
  intrinsic_acceleration_ =
      std::make_unique<IntrinsicAcceleration>(acceleration);
}

template<typename Frame>
void Trajectory<Frame>::clear_intrinsic_acceleration() {
  intrinsic_acceleration_.reset();
}

template<typename Frame>
bool Trajectory<Frame>::has_intrinsic_acceleration() const {
  return intrinsic_acceleration_ != nullptr;
}

template<typename Frame>
Vector<Acceleration, Frame> Trajectory<Frame>::evaluate_intrinsic_acceleration(
    Instant const& time) const {
  if (intrinsic_acceleration_ != nullptr &&
      (fork_ == nullptr || time > fork_->timeline->first)) {
    return (*intrinsic_acceleration_)(time);
  } else {
    return Vector<Acceleration, Frame>({0 * SIUnit<Acceleration>(),
                                        0 * SIUnit<Acceleration>(),
                                        0 * SIUnit<Acceleration>()});
  }
}

template<typename Frame>
void Trajectory<Frame>::WriteToMessage(
    not_null<serialization::Trajectory*> const message) const {
  CHECK(is_root());
  WriteSubTreeToMessage(message);
}

template<typename Frame>
std::unique_ptr<Trajectory<Frame>> Trajectory<Frame>::ReadFromMessage(
    serialization::Trajectory const& message,
    not_null<Body const*> const body) {
  auto trajectory = std::make_unique<Trajectory>(body);
  trajectory->FillSubTreeFromMessage(message);
  return trajectory;
}

template<typename Frame>
void Trajectory<Frame>::WritePointerToMessage(
    not_null<serialization::Trajectory::Pointer*> const message) const {
  not_null<Trajectory const*> ancestor = this;
  while (ancestor->parent_ != nullptr) {
    Fork const& fork = *ancestor->fork_;
    ancestor = ancestor->parent_;
    int const children_distance =
        std::distance(ancestor->children_.begin(), fork.children);
    int const timeline_distance =
        std::distance(ancestor->timeline_.begin(), fork.timeline);
    auto* const fork_message = message->add_fork();
    fork_message->set_children_distance(children_distance);
    fork_message->set_timeline_distance(timeline_distance);
  }
}

template<typename Frame>
not_null<Trajectory<Frame>*> Trajectory<Frame>::ReadPointerFromMessage(
    serialization::Trajectory::Pointer const& message,
    not_null<Trajectory*> const trajectory) {
  CHECK(trajectory->is_root());
  not_null<Trajectory*> descendant = trajectory;
  for (int i = 0; i < message.fork_size(); ++i) {
    auto const& fork_message = message.fork(i);
    int const children_distance = fork_message.children_distance();
    int const timeline_distance = fork_message.timeline_distance();
    auto children_it = descendant->children_.begin();
    auto timeline_it = descendant->timeline_.begin();
    std::advance(children_it, children_distance);
    std::advance(timeline_it, timeline_distance);
    descendant = &children_it->second;
  }
  return descendant;
}

template<typename Frame>
typename Trajectory<Frame>::Iterator&
Trajectory<Frame>::Iterator::operator++() {
  if (!forks_.empty() && current_ == forks_.front().timeline) {
    // Skip over any empty timeline.
    do {
      ancestry_.pop_front();
      forks_.pop_front();
      current_ = ancestry_.front()->timeline_.begin();
    } while (!forks_.empty() && current_ == ancestry_.front()->timeline_.end());
  } else {
    CHECK(current_ != ancestry_.front()->timeline_.end())
        << "Incrementing beyond end of trajectory";
    ++current_;
  }
  CHECK(!current_is_misplaced());
  return *this;
}

template<typename Frame>
bool Trajectory<Frame>::Iterator::at_end() const {
  return forks_.empty() && current_ == ancestry_.front()->timeline_.end();
}

template<typename Frame>
Instant const& Trajectory<Frame>::Iterator::time() const {
  return current_->first;
}

template<typename Frame>
void Trajectory<Frame>::Iterator::InitializeFirst(
    not_null<Trajectory const*> const trajectory) {
  not_null<Trajectory const*> ancestor = trajectory;
  while (ancestor->parent_ != nullptr) {
    ancestry_.push_front(ancestor);
    forks_.push_front(*ancestor->fork_);
    ancestor = ancestor->parent_;
  }
  ancestry_.push_front(ancestor);
  current_ = ancestor->timeline_.begin();
  CHECK(!current_is_misplaced());
}

template<typename Frame>
void Trajectory<Frame>::Iterator::InitializeOnOrAfter(
  Instant const& time, not_null<Trajectory const*> const trajectory) {
  not_null<Trajectory const*> ancestor = trajectory;
  while (ancestor->fork_ != nullptr &&
         time <= ancestor->fork_->timeline->first) {
    ancestry_.push_front(ancestor);
    forks_.push_front(*ancestor->fork_);
    ancestor = ancestor->parent_;
  }
  ancestry_.push_front(ancestor);
  current_ = ancestor->timeline_.lower_bound(time);
  CHECK(!current_is_misplaced());
}

template<typename Frame>
void Trajectory<Frame>::Iterator::InitializeLast(
    not_null<Trajectory const*> const trajectory) {
  // We don't need to keep track of the forks or of the ancestry.
  bool has_empty_timeline = false;
  not_null<Trajectory const*> ancestor = trajectory;
  while (ancestor->timeline_.empty()) {
    has_empty_timeline = true;
    CHECK(ancestor->fork_ != nullptr) << "Empty trajectory";
    current_ = ancestor->fork_->timeline;
    ancestor = ancestor->parent_;
  }
  ancestry_.push_front(ancestor);
  if (!has_empty_timeline) {
    current_ = --ancestor->timeline_.end();
  }
  CHECK(!current_is_misplaced());
}

template<typename Frame>
typename Trajectory<Frame>::Timeline::const_iterator
Trajectory<Frame>::Iterator::current() const {
  return current_;
}

template<typename Frame>
not_null<Trajectory<Frame> const*>
Trajectory<Frame>::Iterator::trajectory() const {
  return ancestry_.back();
}

template<typename Frame>
bool Trajectory<Frame>::Iterator::current_is_misplaced() const {
  return !forks_.empty() && current_ == ancestry_.front()->timeline_.end();
}

template<typename Frame>
DegreesOfFreedom<Frame> const&
Trajectory<Frame>::NativeIterator::degrees_of_freedom() const {
  return this->current()->second;
}

template<typename Frame>
template<typename ToFrame>
DegreesOfFreedom<ToFrame>
Trajectory<Frame>::TransformingIterator<ToFrame>::degrees_of_freedom() const {
  auto it = this->current();
  return transform_(it->first, it->second, this->trajectory());
}

template<typename Frame>
template<typename ToFrame>
Trajectory<Frame>::TransformingIterator<ToFrame>::TransformingIterator(
    Transform<ToFrame> const& transform)
    : Iterator(),
      transform_(transform) {}

template<typename Frame>
Trajectory<Frame>::Trajectory(not_null<Body const*> const body,
                              not_null<Trajectory*> const parent,
                              Fork const& fork)
    : body_(body),
      fork_(new Fork(fork)),
      parent_(parent) {}

template<typename Frame>
bool Trajectory<Frame>::InTimelineOrAtFork(Instant const time) const {
  // Check that |time| is the time of one of our Timeline or the time of fork.
  auto const it = timeline_.find(time);
  if (it == timeline_.end()) {
    if (fork_ == nullptr) {
      return false;
    }

    // Skip over empty timelines to check that this is the fork time.
    Trajectory const* ancestor = parent_;
    Fork fork = *fork_;
    while (ancestor != nullptr && fork.timeline == ancestor->timeline_.end()) {
      fork = *ancestor->fork_;
      ancestor = ancestor->parent_;
    }

    return ancestor != nullptr && fork.timeline->first == time;
  } else {
    return true;
  }
}

template<typename Frame>
void Trajectory<Frame>::WriteSubTreeToMessage(
    not_null<serialization::Trajectory*> const message) const {
  Instant last_instant;
  bool is_first = true;
  serialization::Trajectory::Litter* litter = nullptr;
  for (auto const& pair : children_) {
    Instant const& fork_time = pair.first;
    Trajectory const& child = pair.second;
    if (is_first || fork_time != last_instant) {
      is_first = false;
      last_instant = fork_time;
      litter = message->add_children();
      fork_time.WriteToMessage(litter->mutable_fork_time());
    }
    child.WriteSubTreeToMessage(litter->add_trajectories());
  }
  for (auto const& pair : timeline_) {
    Instant const& instant = pair.first;
    DegreesOfFreedom<Frame> const& degrees_of_freedom = pair.second;
    auto const instantaneous_degrees_of_freedom = message->add_timeline();
    instant.WriteToMessage(instantaneous_degrees_of_freedom->mutable_instant());
    degrees_of_freedom.WriteToMessage(
        instantaneous_degrees_of_freedom->mutable_degrees_of_freedom());
  }
}

template<typename Frame>
void Trajectory<Frame>::FillSubTreeFromMessage(
    serialization::Trajectory const& message) {
  auto timeline_it = message.timeline().begin();
  for (serialization::Trajectory::Litter const& litter : message.children()) {
    Instant const fork_time = Instant::ReadFromMessage(litter.fork_time());
    for (;
         timeline_it != message.timeline().end() &&
         Instant::ReadFromMessage(timeline_it->instant()) <= fork_time;
         ++timeline_it) {
      Append(Instant::ReadFromMessage(timeline_it->instant()),
             DegreesOfFreedom<Frame>::ReadFromMessage(
                 timeline_it->degrees_of_freedom()));
    }
    for (serialization::Trajectory const& child : litter.trajectories()) {
      NewFork(fork_time)->FillSubTreeFromMessage(child);
    }
  }
  for (; timeline_it != message.timeline().end(); ++timeline_it) {
    Append(Instant::ReadFromMessage(timeline_it->instant()),
           DegreesOfFreedom<Frame>::ReadFromMessage(
               timeline_it->degrees_of_freedom()));
  }
}

}  // namespace physics
}  // namespace principia
