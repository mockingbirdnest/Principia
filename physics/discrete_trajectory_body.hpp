#pragma once

#include "physics/discrete_trajectory.hpp"

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
DiscreteTrajectory<Frame>::DiscreteTrajectory(not_null<Body const*> const body)
    : body_(body) {
  CHECK(body_->is_compatible_with<Frame>())
      << "Oblate body not in the same frame as the trajectory";
}

template<typename Frame>
DiscreteTrajectory<Frame>::~DiscreteTrajectory() {
  if (on_destroy_) {
    on_destroy_(this);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::set_on_destroy(
    std::function<void(not_null<DiscreteTrajectory<Frame>const *> const)>
        on_destroy) {
  on_destroy_ = on_destroy;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::NativeIterator
DiscreteTrajectory<Frame>::first() const {
  NativeIterator it = Begin();
  return it;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::NativeIterator
DiscreteTrajectory<Frame>::on_or_after(
    Instant const& time) const {
  NativeIterator it = Find(time);//TODO(phl):Not quite the same.
  return it;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::NativeIterator
DiscreteTrajectory<Frame>::last() const {
  NativeIterator it = End();
  return --it;
}

template<typename Frame>
template<typename ToFrame>
typename DiscreteTrajectory<Frame>::TEMPLATE TransformingIterator<ToFrame>
DiscreteTrajectory<Frame>::first_with_transform(
    Transform<ToFrame> const& transform) const {
  TransformingIterator<ToFrame> it(transform);
  it.InitializeFirst(this);
  return it;
}

template<typename Frame>
template<typename ToFrame>
typename DiscreteTrajectory<Frame>::TEMPLATE TransformingIterator<ToFrame>
DiscreteTrajectory<Frame>::on_or_after_with_transform(
    Instant const& time,
    Transform<ToFrame> const& transform) const {
  TransformingIterator<ToFrame> it(transform);
  it.InitializeOnOrAfter(time, this);
  return it;
}

template<typename Frame>
template<typename ToFrame>
typename DiscreteTrajectory<Frame>::TEMPLATE TransformingIterator<ToFrame>
DiscreteTrajectory<Frame>::last_with_transform(
    Transform<ToFrame> const& transform) const {
  TransformingIterator<ToFrame> it(transform);
  it.InitializeLast(this);
  return it;
}

template<typename Frame>
std::map<Instant, Position<Frame>>
DiscreteTrajectory<Frame>::Positions() const {
  std::map<Instant, Position<Frame>> result;
  for (auto it = Begin(); it != End(); ++it) {
    auto timeline_it = it.current();
    Instant const& time = timeline_it->first;
    DegreesOfFreedom<Frame> const& degrees_of_freedom = timeline_it->second;
    result.emplace_hint(result.end(), time, degrees_of_freedom.position());
  }
  return result;
}

template<typename Frame>
std::map<Instant, Velocity<Frame>> DiscreteTrajectory<Frame>::Velocities() const {
  std::map<Instant, Velocity<Frame>> result;
  for (auto it = Begin(); it != End(); ++it) {
    auto timeline_it = it.current();
    Instant const& time = timeline_it->first;
    DegreesOfFreedom<Frame> const& degrees_of_freedom = timeline_it->second;
    result.emplace_hint(result.end(), time, degrees_of_freedom.velocity());
  }
  return result;
}

template<typename Frame>
std::list<Instant> DiscreteTrajectory<Frame>::Times() const {
  std::list<Instant> result;
  for (auto it = Begin(); it != End(); ++it) {
    auto timeline_it = it.current();
    Instant const& time = timeline_it->first;
    result.push_back(time);
  }
  return result;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  if (!first().at_end() && last().time() == time) {
    LOG(WARNING) << "Append at existing time " << time
                 << ", time range = [" << Times().front() << ", "
                 << Times().back() << "]";
    return;
  }
  auto it = timeline_.emplace_hint(timeline_.end(),
                                   time,
                                   degrees_of_freedom);
  CHECK(timeline_.end() == ++it) << "Append out of order";
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetAfter(Instant const& time) {
  // Each of these blocks gets an iterator denoting the first entry with
  // time > |time|.  It then removes that entry and all the entries that follow
  // it.  This preserve any entry with time == |time|.
  {
    auto const it = timeline_.upper_bound(time);
    CHECK(is_root() || time >= ForkTime())
        << "ForgetAfter before the fork time";
    timeline_.erase(it, timeline_.end());
  }
  {
    auto const it = children_.upper_bound(time);
    children_.erase(it, children_.end());
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetBefore(Instant const& time) {
  // Check that this is a root.
  CHECK(is_root()) << "ForgetBefore on a nonroot trajectory";
  // Each of these blocks gets an iterator denoting the first entry with
  // time > |time|.  It then removes that all the entries that precede it.  This
  // removes any entry with time == |time|.
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
template<typename B>
std::enable_if_t<std::is_base_of<Body, B>::value, not_null<B const*>>
DiscreteTrajectory<Frame>::body() const {
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
void DiscreteTrajectory<Frame>::set_intrinsic_acceleration(
    IntrinsicAcceleration const acceleration) {
  CHECK(body_->is_massless()) << "DiscreteTrajectory is for a massive body";
  CHECK(intrinsic_acceleration_ == nullptr)
      << "DiscreteTrajectory already has an intrinsic acceleration";
  intrinsic_acceleration_ =
      std::make_unique<IntrinsicAcceleration>(acceleration);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::clear_intrinsic_acceleration() {
  intrinsic_acceleration_.reset();
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::has_intrinsic_acceleration() const {
  return intrinsic_acceleration_ != nullptr;
}

template<typename Frame>
Vector<Acceleration, Frame>
DiscreteTrajectory<Frame>::evaluate_intrinsic_acceleration(
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
void DiscreteTrajectory<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> const message) const {
  LOG(INFO) << __FUNCTION__;
  CHECK(is_root());
  WriteSubTreeToMessage(message);
  LOG(INFO) << NAMED(this);
  LOG(INFO) << NAMED(message->SpaceUsed());
  LOG(INFO) << NAMED(message->ByteSize());
}

template<typename Frame>
std::unique_ptr<DiscreteTrajectory<Frame>>
DiscreteTrajectory<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectory const& message,
    not_null<Body const*> const body) {
  auto trajectory = std::make_unique<DiscreteTrajectory>(body);
  trajectory->FillSubTreeFromMessage(message);
  return trajectory;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WritePointerToMessage(
    not_null<serialization::DiscreteTrajectory::Pointer*> const message) const {
  not_null<DiscreteTrajectory const*> ancestor = this;
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
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::ReadPointerFromMessage(
    serialization::DiscreteTrajectory::Pointer const& message,
    not_null<DiscreteTrajectory*> const trajectory) {
  CHECK(trajectory->is_root());
  not_null<DiscreteTrajectory*> descendant = trajectory;
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
typename DiscreteTrajectory<Frame>::Iterator&
DiscreteTrajectory<Frame>::Iterator::operator++() {
  if (!forks_.empty() && current_ == forks_.front().timeline) {
    // Skip over any timeline where the fork is at |end()|.  These are the ones
    // that were forked at the fork point of their parent.  Looking at the
    // |begin()| of the parent would be wrong (the fork would see changes to its
    // parent after the fork point).
    do {
      ancestry_.pop_front();
      forks_.pop_front();
    } while (!forks_.empty() &&
             forks_.front().timeline == ancestry_.front()->timeline_.end());
    current_ = ancestry_.front()->timeline_.begin();
  } else {
    CHECK(current_ != ancestry_.front()->timeline_.end())
        << "Incrementing beyond end of trajectory";
    ++current_;
  }
  CHECK(!current_is_misplaced());
  return *this;
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::Iterator::at_end() const {
  return forks_.empty() && current_ == ancestry_.front()->timeline_.end();
}

template<typename Frame>
Instant const& DiscreteTrajectory<Frame>::Iterator::time() const {
  return current_->first;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame> const*>
DiscreteTrajectory<Frame>::Iterator::trajectory() const {
  return ancestry_.back();
}

template<typename Frame>
DegreesOfFreedom<Frame> const&
DiscreteTrajectory<Frame>::NativeIterator::degrees_of_freedom() const {
  return this->current()->second;
}

template<typename Frame>
template<typename ToFrame>
DegreesOfFreedom<ToFrame>
DiscreteTrajectory<Frame>::
TransformingIterator<ToFrame>::degrees_of_freedom() const {
  auto it = this->current();
  return transform_(it->first, it->second, this->trajectory());
}

template<typename Frame>
template<typename ToFrame>
DiscreteTrajectory<Frame>::TransformingIterator<ToFrame>::TransformingIterator(
    Transform<ToFrame> const& transform)
    : Iterator(),
      transform_(transform) {}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*> DiscreteTrajectory<Frame>::that() {
  return this;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame> const*>
DiscreteTrajectory<Frame>::that() const {
  return this;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_begin() const {
  return timeline_.begin();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_end() const {
  return timeline_.end();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_find(Instant const& time) const {
  return timeline_.find(time);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::timeline_insert(TimelineConstIterator begin,
                                                TimelineConstIterator end) {
  return timeline_.insert(begin, end);
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::timeline_empty() const {
  return timeline_.empty();
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteSubTreeToMessage(
    not_null<serialization::DiscreteTrajectory*> const message) const {
  Instant last_instant;
  bool is_first = true;
  serialization::DiscreteTrajectory::Litter* litter = nullptr;
  for (auto const& pair : children_) {
    Instant const& fork_time = pair.first;
    DiscreteTrajectory const& child = pair.second;
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
void DiscreteTrajectory<Frame>::FillSubTreeFromMessage(
    serialization::DiscreteTrajectory const& message) {
  auto timeline_it = message.timeline().begin();
  for (serialization::DiscreteTrajectory::Litter const& litter : message.children()) {
    Instant const fork_time = Instant::ReadFromMessage(litter.fork_time());
    for (;
         timeline_it != message.timeline().end() &&
         Instant::ReadFromMessage(timeline_it->instant()) <= fork_time;
         ++timeline_it) {
      Append(Instant::ReadFromMessage(timeline_it->instant()),
             DegreesOfFreedom<Frame>::ReadFromMessage(
                 timeline_it->degrees_of_freedom()));
    }
    for (serialization::DiscreteTrajectory const& child : litter.trajectories()) {
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
