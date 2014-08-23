#pragma once

#include "trajectory.hpp"

#include <list>
#include <map>

#ifndef _MANAGED
#include "glog/logging.h"
#endif

namespace principia {
namespace physics {

template<typename Frame>
Trajectory<Frame>::Trajectory(Body const* body)
    :
#ifdef _MANAGED
      body_(body),
#else
      body_(CHECK_NOTNULL(body)),
#endif
      parent_(nullptr) {}

template<typename Frame>
std::map<Time, Vector<Length, Frame>> Trajectory<Frame>::Positions() const {
  return ApplyToDegreesOfFreedom<Vector<Length, Frame>>(
      [](DegreesOfFreedom<Frame> const& s) { return s.position; });
}

template<typename Frame>
std::map<Time, Vector<Speed, Frame>> Trajectory<Frame>::Velocities() const {
  return ApplyToDegreesOfFreedom<Vector<Speed, Frame>>(
      [](DegreesOfFreedom<Frame> const& s) { return s.velocity; });
}

template<typename Frame>
std::list<Time> Trajectory<Frame>::Times() const {
  std::list<Time> result;

  // Our own data points in increasing time order.
  for (const auto it : timeline_) {
    Time const& time = it.first;
    result.push_back(time);
  }

  // The data points of our ancestors in decreasing time order.
  Trajectory const* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    Timeline::iterator it = *ancestor->fork_;
    do {
      Time const& time = it->first;
      result.push_front(time);
    } while (it-- !=  // Postdecrement to process begin.
             ancestor->parent_->timeline_.begin());
    ancestor = ancestor->parent_;
  }

  return result;
}

template<typename Frame>
Vector<Length, Frame> const& Trajectory<Frame>::last_position() const {
  return timeline_.rbegin()->second.position;
}

template<typename Frame>
Vector<Speed, Frame> const& Trajectory<Frame>::last_velocity() const {
  return timeline_.rbegin()->second.velocity;
}

template<typename Frame>
Time const& Trajectory<Frame>::last_time() const {
  return timeline_.rbegin()->first;
}

template<typename Frame>
void Trajectory<Frame>::Append(Vector<Length, Frame> const& position,
                               Vector<Speed, Frame> const& velocity,
                               Time const& time) {
  auto inserted = timeline_.insert(
      std::make_pair(time, DegreesOfFreedom<Frame>(position, velocity)));
#ifndef _MANAGED
  CHECK(timeline_.end() == ++inserted.first) << "Append out of order";
  CHECK(inserted.second) << "Append at existing time";
#endif
}

template<typename Frame>
void Trajectory<Frame>::ForgetAfter(Time const& time) {
  // Check that |time| is the time of one of our Timeline or the time of fork.
  auto const it = timeline_.find(time);
  if (it == timeline_.end()) {
#ifndef _MANAGED
    CHECK(fork_ != nullptr)
        << "ForgetAfter a nonexistent time for a root trajectory";
    CHECK_EQ((*fork_)->first, time)
        << "ForgetAfter a nonexistent time for a nonroot trajectory";
#endif
  }

  // Each of these blocks gets an iterator denoting the first entry with
  // time > |time|.  It then removes that entry and all the entries that follow
  // it.  This preserve any entry with time == |time|.
  {
    const auto it = timeline_.upper_bound(time);
    timeline_.erase(it, timeline_.end());
  }
  {
    const auto it = children_.upper_bound(time);
    children_.erase(it, children_.end());
  }
}

template<typename Frame>
void Trajectory<Frame>::ForgetBefore(Time const& time) {
#ifndef _MANAGED
  // Check that this is a root.
  CHECK(is_root()) << "ForgetBefore on a nonroot trajectory";
  // Check that |time| is the time of one of our Timeline or the time of fork.
  CHECK(timeline_.find(time) != timeline_.end())
      << "ForgetBefore a nonexistent time";
#endif
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
Trajectory<Frame>* Trajectory<Frame>::Fork(Time const& time) {
  auto fork_it = timeline_.find(time);
#ifndef _MANAGED
  CHECK(fork_it != timeline_.end()) << "Fork at nonexistent time";
#endif
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
Time const* Trajectory<Frame>::fork_time() const {
  if (parent_ == nullptr) {
    return nullptr;
  } else {
    return &((*fork_)->first);
  }
}

template<typename Frame>
Body const& Trajectory<Frame>::body() const {
  return *body_;
}

template<typename Frame>
Trajectory<Frame>::Trajectory(Body const* const body,
                              Trajectory* const parent,
                              typename Timeline::iterator const& fork)
    :
#ifdef _MANAGED
      body_(body),
      parent_(parent),
#else
      body_(CHECK_NOTNULL(body)),
      parent_(CHECK_NOTNULL(parent)),
#endif
      fork_(new Timeline::iterator(fork)) {}

template<typename Frame>
template<typename Value>
std::map<Time, Value> Trajectory<Frame>::ApplyToDegreesOfFreedom(
    std::function<Value(DegreesOfFreedom<Frame> const&)> compute_value) const {
  std::map<Time, Value> result;

  // Our own data points in increasing time order.
  for (const auto it : timeline_) {
    Time const& time = it.first;
    DegreesOfFreedom<Frame> const& degrees_of_freedom = it.second;
    result.insert(result.end(),
                  std::make_pair(time, compute_value(degrees_of_freedom)));
  }

  // The data points of our ancestors in decreasing time order.
  Trajectory const* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    Timeline::iterator it = *ancestor->fork_;
    do {
      Time const& time = it->first;
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
