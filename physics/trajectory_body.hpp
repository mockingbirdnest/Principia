#pragma once

#include "trajectory.hpp"

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
  return GetState<Vector<Length, Frame>>(
      [](State const& s) { return s.position(); });
}

template<typename Frame>
std::map<Time, Vector<Speed, Frame>> Trajectory<Frame>::Velocities() const {
  return GetState<Vector<Speed, Frame>>(
      [](State const& s) { return s.velocity(); });
}

template<typename Frame>
std::list<Time> Trajectory<Frame>::Times() const {
  std::list<Time> result;

  // Our own data points in increasing time order.
  for (const auto it : states_) {
    Time const& time = it.first;
    result.push_back(time);
  }

  // The data points of our ancestors in decreasing time order.
  Trajectory const* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    States::iterator it = *ancestor->parent_state_;
    do {
      Time const& time = it->first;
      result.push_front(time);
    } while (it-- !=  // Postdecrement to process begin.
             ancestor->parent_->states_.begin());
    ancestor = ancestor->parent_;
  }

  return result;
}

template<typename Frame>
void Trajectory<Frame>::Append(Vector<Length, Frame> const& position,
                               Vector<Speed, Frame> const& velocity,
                               Time const& time) {
  auto inserted =
      states_.insert(std::make_pair(time, State(position, velocity)));
#ifndef _MANAGED
  CHECK(states_.end() == ++inserted.first) << "Append out of order";
  CHECK(inserted.second) << "Append at existing time";
#endif
}

template<typename Frame>
void Trajectory<Frame>::ForgetAfter(Time const& time) {
  // Each of these blocks gets an iterator denoting the first entry with
  // time > |time|.  It then removes that entry and all the entries that follow
  // it.  This preserve any entry with time == |time|.
  {
    const auto it = states_.upper_bound(time);
    states_.erase(it, states_.end());
  }
  {
    const auto it = children_.upper_bound(time);
    children_.erase(it, children_.end());
  }
  {
    const auto it = bursts_.upper_bound(time);
    burst_.erase(it, bursts_.end());
  }
}

template<typename Frame>
void Trajectory<Frame>::ForgetBefore(Time const& time) {
  {
    const auto it = states_.lower_bound(time);
    states_.erase(states_.begin(), it);
  }
  {
    const auto it = children_.lower_bound(time);
    children_.erase(children_.begin(), it);
  }
  //TODO(phl): This is not correct because there may be a burst which starts
  // before |time| but ends after.
  {
    const auto it = bursts_.lower_bound(time);
    burst_.erase(bursts_.begin(), it);
  }
}

template<typename Frame>
Trajectory<Frame>* Trajectory<Frame>::Fork(Time const& time) {
  auto fork_it = states_.find(time);
  CHECK(fork_it != states_.end()) << "Fork at nonexistent time";
  std::unique_ptr<Trajectory<Frame>> child(
      new Trajectory(body_, this /*parent*/, fork_it));
  child->states_.insert(++fork_it, states_.end());
  auto const child_it = children_.emplace(time, std::move(child));
  return child_it->second.get();
}

template<typename Frame>
Body const* Trajectory<Frame>::body() const {
  return body_;
}

template<typename Frame>
void Trajectory<Frame>::AddBurst(
    Vector<Acceleration, Frame> const& acceleration,
    Time const& time1,
    Time const& time2) {
  Time const duration = time2 - time1;
  if (duration <= 0) {
    return;
  }

  // Check that the times don't straddle an existing point.  This needs some
  // serious testing as I am too tired to be sure that this test is correct.
  const auto it1 = bursts_.upper_bound(time1);
  const auto it2 = bursts_.lower_bound(time2);
  CHECK_EQ(it1, it2);

  bursts_.insert({time1, Burst(acceleration, duration)});
}

template<typename Frame>
Trajectory<Frame>::Trajectory(Body const* const body,
                              Trajectory const* const parent,
                              typename States::iterator const& parent_state)
    : body_(body),
      parent_(parent),//TODO(phl):CHECK_NOTNULL
      parent_state_(new States::iterator(parent_state)) {};

template<typename Frame>
Trajectory<Frame>::Burst::Burst(Vector<Acceleration, Frame> const& acceleration,
                                Time const& duration)
    : acceleration_(acceleration),
      duration_(duration) {}

template<typename Frame>
Trajectory<Frame>::State::State(Vector<Length, Frame> const& position,
                                Vector<Speed, Frame> const& velocity)
    : position_(position),
      velocity_(velocity) {}

template<typename Frame>
Vector<Length, Frame> const& Trajectory<Frame>::State::position() const {
  return position_;
}

template<typename Frame>
Vector<Speed, Frame> const& Trajectory<Frame>::State::velocity() const {
  return velocity_;
}

template<typename Frame>
template<typename Value>
std::map<Time, Value> Trajectory<Frame>::GetState(
    std::function<Value(State const&)> fun) const {
  std::map<Time, Value> result;

  // Our own data points in increasing time order.
  for (const auto it : states_) {
    Time const& time = it.first;
    State const& state = it.second;
    result.insert(result.end(), std::make_pair(time, fun(state)));
  }

  // The data points of our ancestors in decreasing time order.
  Trajectory const* ancestor = this;
  while (ancestor->parent_ != nullptr) {
    States::iterator it = *ancestor->parent_state_;
    do {
      Time const& time = it->first;
      State const& state = it->second;
      result.insert(result.begin(), std::make_pair(time, fun(state)));
    } while (it-- !=  // Postdecrement to process begin.
             ancestor->parent_->states_.begin());
    ancestor = ancestor->parent_;
  }

  return result;
}

}  // namespace physics
}  // namespace principia
