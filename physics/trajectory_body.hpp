#pragma once

#include "trajectory.hpp"

#include <list>
#include <map>

#include "geometry/named_quantities.hpp"
#include "glog/logging.h"

using principia::geometry::Instant;

namespace principia {
namespace physics {

template<typename Frame>
Trajectory<Frame>::Trajectory(Body<Frame> const& body)
    : body_(body),
      parent_(nullptr) {}

template<typename Frame>
typename Trajectory<Frame>::NativeIterator Trajectory<Frame>::first() const {
  NativeIterator it;
  it.InitializeFirst(this);
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
typename Trajectory<Frame>::TransformingIterator<ToFrame>
Trajectory<Frame>::first_with_transform(
    Transform<ToFrame> const& transform) const {
  TransformingIterator<ToFrame> it(transform);
  it.InitializeFirst(this);
  return it;
}

template<typename Frame>
template<typename ToFrame>
typename Trajectory<Frame>::TransformingIterator<ToFrame>
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
    result.insert(result.end(),
                  std::make_pair(time, it.degrees_of_freedom().position));
  }
  return result;
}

template<typename Frame>
std::map<Instant, Velocity<Frame>> Trajectory<Frame>::Velocities() const {
  std::map<Instant, Velocity<Frame>> result;
  for (NativeIterator it = first(); !it.at_end(); ++it) {
    Instant const& time = it.time();
    result.insert(result.end(),
                  std::make_pair(time, it.degrees_of_freedom().velocity));
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
const DegreesOfFreedom<Frame>& Trajectory<Frame>::GetDegreesOfFreedom(
    Instant const& time) const {
  Trajectory const* ancestor = this;
  while (ancestor->fork_ != nullptr && time <= (*ancestor->fork_)->first) {
    ancestor = ancestor->parent_;
  }
  auto const it = ancestor->timeline_.find(time);
  CHECK(it != ancestor->timeline_.end())
      << "Time " << time << " not in trajectory";
  return it->second;
}

template<typename Frame>
void Trajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  auto inserted = timeline_.emplace(time, degrees_of_freedom);
  CHECK(timeline_.end() == ++inserted.first) << "Append out of order";
  CHECK(inserted.second) << "Append at existing time";
}

template<typename Frame>
void Trajectory<Frame>::ForgetAfter(Instant const& time) {
  // Check that |time| is the time of one of our Timeline or the time of fork.
  auto const it = timeline_.find(time);
  if (it == timeline_.end()) {
    CHECK(fork_ != nullptr)
        << "ForgetAfter a nonexistent time for a root trajectory";
    CHECK_EQ((*fork_)->first, time)
        << "ForgetAfter a nonexistent time for a nonroot trajectory";
  }

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
  // Check that |time| is the time of one of our Timeline or the time of fork.
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
Trajectory<Frame>* Trajectory<Frame>::Fork(Instant const& time) {
  auto fork_it = timeline_.find(time);
  CHECK(fork_it != timeline_.end()) << "Fork at nonexistent time";
  std::unique_ptr<Trajectory<Frame>> child(
      new Trajectory(body_, this /*parent*/, fork_it));
  child->timeline_.insert(++fork_it, timeline_.end());
  auto const child_it = children_.emplace(time, std::move(child));
  return child_it->second.get();
}

template<typename Frame>
void Trajectory<Frame>::DeleteFork(Trajectory** const fork) {
  CHECK_NOTNULL(fork);
  CHECK_NOTNULL(*fork);
  Instant const* const fork_time = (*fork)->fork_time();
  CHECK_NOTNULL(fork_time);
  // Find the position of |*fork| among our children and remove it.
  auto const range = children_.equal_range(*fork_time);
  for (auto it = range.first; it != range.second; ++it) {
    if (it->second.get() == *fork) {
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
Trajectory<Frame> const* Trajectory<Frame>::root() const {
  Trajectory const* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    ancestor = ancestor->parent_;
  }
  return ancestor;
}

template<typename Frame>
Trajectory<Frame>* Trajectory<Frame>::root() {
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
    return &((*fork_)->first);
  }
}

template<typename Frame>
Body<Frame> const& Trajectory<Frame>::body() const {
  return body_;
}

template<typename Frame>
void Trajectory<Frame>::set_intrinsic_acceleration(
    IntrinsicAcceleration const acceleration) {
  CHECK(body_.is_massless()) << "Trajectory is for a massive body";
  CHECK(intrinsic_acceleration_ == nullptr)
      << "Trajectory already has an intrinsic acceleration";
  intrinsic_acceleration_.reset(new IntrinsicAcceleration(acceleration));
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
      (fork_ == nullptr || time > (*fork_)->first)) {
    return (*intrinsic_acceleration_)(time);
  } else {
    return Vector<Acceleration, Frame>({0 * SIUnit<Acceleration>(),
                                        0 * SIUnit<Acceleration>(),
                                        0 * SIUnit<Acceleration>()});
  }
}

template<typename Frame>
typename Trajectory<Frame>::Iterator&
Trajectory<Frame>::Iterator::operator++() {
  if (!forks_.empty() && current_ == forks_.front()) {
    ancestry_.pop_front();
    forks_.pop_front();
    current_ = ancestry_.front()->timeline_.begin();
  } else {
    CHECK(current_ != ancestry_.front()->timeline_.end());
    ++current_;
  }
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
    Trajectory const* trajectory) {
  CHECK_NOTNULL(trajectory);
  Trajectory const* ancestor = trajectory;
  while (ancestor->parent_ != nullptr) {
    ancestry_.push_front(ancestor);
    forks_.push_front(*ancestor->fork_);
    ancestor = ancestor->parent_;
  }
  ancestry_.push_front(ancestor);
  current_ = ancestor->timeline_.begin();
}

template<typename Frame>
void Trajectory<Frame>::Iterator::InitializeLast(
    Trajectory const* trajectory) {
  CHECK_NOTNULL(trajectory);
  // We don't need to really keep track of the forks or of the ancestry.
  if (trajectory->timeline_.empty()) {
    CHECK(trajectory->fork_ != nullptr) << "Empty trajectory";
    ancestry_.push_front(trajectory->parent_);
    current_ = *trajectory->fork_;
  } else {
    ancestry_.push_front(trajectory);
    current_ = --trajectory->timeline_.end();
  }
}

template<typename Frame>
typename Trajectory<Frame>::Timeline::const_iterator
Trajectory<Frame>::Iterator::current() const {
  return current_;
}

template<typename Frame>
DegreesOfFreedom<Frame> const&
Trajectory<Frame>::NativeIterator::degrees_of_freedom() const {
  return this->current()->second;
}

template<typename Frame>
template<typename ToFrame>
DegreesOfFreedom<ToFrame> const
Trajectory<Frame>::TransformingIterator<ToFrame>::degrees_of_freedom() const {
  auto it = current();
  return transform_(it->first, it->second);
}

template<typename Frame>
template<typename ToFrame>
Trajectory<Frame>::TransformingIterator<ToFrame>::TransformingIterator(
    Transform<ToFrame> const& transform)
    : Iterator(),
      transform_(transform) {}

template<typename Frame>
Trajectory<Frame>::Trajectory(Body<Frame> const& body,
                              Trajectory* const parent,
                              typename Timeline::iterator const& fork)
    : body_(body),
      parent_(CHECK_NOTNULL(parent)),
      fork_(new typename Timeline::iterator(fork)) {}

}  // namespace physics
}  // namespace principia
