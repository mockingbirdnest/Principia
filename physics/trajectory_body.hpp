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
Trajectory<Frame>::Trajectory(Body const& body)
    : body_(body),
      parent_(nullptr) {}

template<typename Frame>
Trajectory<Frame>::~Trajectory() {
  Instant const* const fork_time = fork_time();
  if (fork_time != nullptr) {
    // Find our position among our siblings and remove ourselves.
    auto siblings = parent_->children_;
    auto const range = siblings->equal_range(*fork_time);
    for (auto it = range.first; it != range.second; ++it) {
      if (it->second.get() == this) {
        // I now own myself for a few nanoseconds.
        it->second.release();
        parent_->children_->erase(it);
        return;
      }
    }
    LOG(FATAL) << "Inconsistent parent/chidren";
  }
}

template<typename Frame>
std::map<Instant, Position<Frame>>
Trajectory<Frame>::Positions() const {
  return ApplyToDegreesOfFreedom<Position<Frame>>(
      [](DegreesOfFreedom<Frame> const& s) { return s.position; });
}

template<typename Frame>
std::map<Instant, Velocity<Frame>> Trajectory<Frame>::Velocities() const {
  return ApplyToDegreesOfFreedom<Velocity<Frame>>(
      [](DegreesOfFreedom<Frame> const& s) { return s.velocity; });
}

template<typename Frame>
std::list<Instant> Trajectory<Frame>::Times() const {
  std::list<Instant> result;

  // Our own data points in increasing time order.
  for (auto const it : timeline_) {
    Instant const& time = it.first;
    result.push_back(time);
  }

  // The data points of our ancestors in decreasing time order.
  Trajectory const* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    Timeline::iterator it = *ancestor->fork_;
    do {
      Instant const& time = it->first;
      result.push_front(time);
    } while (it-- !=  // Postdecrement to process begin.
             ancestor->parent_->timeline_.begin());
    ancestor = ancestor->parent_;
  }

  return result;
}

template<typename Frame>
Position<Frame> const& Trajectory<Frame>::last_position() const {
  if (timeline_.empty()) {
    CHECK(fork_ != nullptr) << "Empty trajectory";
    return (*fork_)->second.position;
  } else {
    return timeline_.rbegin()->second.position;
  }
}

template<typename Frame>
Velocity<Frame> const& Trajectory<Frame>::last_velocity() const {
  if (timeline_.empty()) {
    CHECK(fork_ != nullptr) << "Empty trajectory";
    return (*fork_)->second.velocity;
  } else {
    return timeline_.rbegin()->second.velocity;
  }
}

template<typename Frame>
Instant const& Trajectory<Frame>::last_time() const {
  if (timeline_.empty()) {
    CHECK(fork_ != nullptr) << "Empty trajectory";
    return (*fork_)->first;
  } else {
    return timeline_.rbegin()->first;
  }
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
Body const& Trajectory<Frame>::body() const {
  return body_;
}

template<typename Frame>
void Trajectory<Frame>::set_intrinsic_acceleration(
    IntrinsicAcceleration&& acceleration) {
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
Trajectory<Frame>::Trajectory(Body const& body,
                              Trajectory* const parent,
                              typename Timeline::iterator const& fork)
    : body_(body),
      parent_(CHECK_NOTNULL(parent)),
      fork_(new Timeline::iterator(fork)) {}

template<typename Frame>
template<typename Value>
std::map<Instant, Value> Trajectory<Frame>::ApplyToDegreesOfFreedom(
    std::function<Value(DegreesOfFreedom<Frame> const&)> compute_value) const {
  std::map<Instant, Value> result;

  // Our own data points in increasing time order.
  for (auto const it : timeline_) {
    Instant const& time = it.first;
    DegreesOfFreedom<Frame> const& degrees_of_freedom = it.second;
    result.insert(result.end(),
                  std::make_pair(time, compute_value(degrees_of_freedom)));
  }

  // The data points of our ancestors in decreasing time order.
  Trajectory const* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    Timeline::iterator it = *ancestor->fork_;
    do {
      Instant const& time = it->first;
      DegreesOfFreedom<Frame> const& degrees_of_freedom = it->second;
      result.insert(result.begin(),
                    std::make_pair(time, compute_value(degrees_of_freedom)));
    } while (it-- !=  // Postdecrement to process begin.
             ancestor->parent_->timeline_.begin());
    ancestor = ancestor->parent_;
  }

  return result;
}

}  // namespace physics
}  // namespace principia
