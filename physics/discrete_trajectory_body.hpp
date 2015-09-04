#pragma once

#include "physics/discrete_trajectory.hpp"

#include <algorithm>
#include <list>
#include <map>

#include "geometry/named_quantities.hpp"
#include "glog/logging.h"

namespace principia {

using base::make_not_null_unique;
using geometry::Instant;

namespace physics {

template<typename Frame>
Instant const& ForkableTraits<DiscreteTrajectory<Frame>>::time(
    TimelineConstIterator const it) {
  return it->first
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
  return NativeIterator(Begin());
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::NativeIterator
DiscreteTrajectory<Frame>::on_or_after(
    Instant const& time) const {
  return NativeIterator(Find(time));//TODO(phl):Not quite the same.
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::NativeIterator
DiscreteTrajectory<Frame>::last() const {
  auto it = End();
  return NativeIterator(--it);
}

template<typename Frame>
template<typename ToFrame>
typename DiscreteTrajectory<Frame>::TEMPLATE TransformingIterator<ToFrame>
DiscreteTrajectory<Frame>::first_with_transform(
    Transform<ToFrame> const& transform) const {
  return TransformingIterator<ToFrame>(Begin(), transform);
}

template<typename Frame>
template<typename ToFrame>
typename DiscreteTrajectory<Frame>::TEMPLATE TransformingIterator<ToFrame>
DiscreteTrajectory<Frame>::on_or_after_with_transform(
    Instant const& time,
    Transform<ToFrame> const& transform) const {
  return TransformingIterator<ToFrame>(Find(time), transform);//TODO(phl):Not quite the same.
}

template<typename Frame>
template<typename ToFrame>
typename DiscreteTrajectory<Frame>::TEMPLATE TransformingIterator<ToFrame>
DiscreteTrajectory<Frame>::last_with_transform(
    Transform<ToFrame> const& transform) const {
  auto it = End();
  return TransformingIterator<ToFrame>(--it, transform);
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
    CHECK(is_root() || time >= *ForkTime())
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
void DiscreteTrajectory<Frame>::set_intrinsic_acceleration(
    IntrinsicAcceleration const acceleration) {
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
      (is_root() || time >= *ForkTime())) {
    return (*intrinsic_acceleration_)(time);
  } else {
    return Vector<Acceleration, Frame>({0 * SIUnit<Acceleration>(),
                                        0 * SIUnit<Acceleration>(),
                                        0 * SIUnit<Acceleration>()});
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteToMessage(
    not_null<serialization::Trajectory*> const message) const {
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
    serialization::Trajectory const& message) {
  auto trajectory = std::make_unique<DiscreteTrajectory>();
  trajectory->FillSubTreeFromMessage(message);
  return trajectory;
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::NativeIterator::at_end() const {
  return *this == trajectory()->End();
}

template<typename Frame>
Instant const& DiscreteTrajectory<Frame>::NativeIterator::time() const {
  return current()->first;
}

template<typename Frame>
typename DegreesOfFreedom<Frame> const&
DiscreteTrajectory<Frame>::NativeIterator::degrees_of_freedom() const {
  return current()->second;
}

template<typename Frame>
DiscreteTrajectory<Frame>::NativeIterator::NativeIterator(Iterator it)
    : Iterator(std::move(it)) {}

template<typename Frame>
template<typename ToFrame>
bool DiscreteTrajectory<Frame>::TransformingIterator<ToFrame>::at_end() const {
  return *this == trajectory()->End();
}

template<typename Frame>
template<typename ToFrame>
Instant const&
DiscreteTrajectory<Frame>::TransformingIterator<ToFrame>::time() const {
  return current()->first;
}
template<typename Frame>
template<typename ToFrame>
typename DegreesOfFreedom<ToFrame>
DiscreteTrajectory<Frame>::
TransformingIterator<ToFrame>::degrees_of_freedom() const {
  auto it = current();
  return transform_(it->first, it->second, trajectory());
}

template<typename Frame>
template<typename ToFrame>
DiscreteTrajectory<Frame>::TransformingIterator<ToFrame>::TransformingIterator(
    Iterator it, Transform<ToFrame> transform)
    : Iterator(std::move(it)),
      transform_(std::move(transform)) {}

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
    not_null<serialization::Trajectory*> const message) const {
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
    serialization::Trajectory const& message) {
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
