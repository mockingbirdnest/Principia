
#pragma once

#include "physics/discrete_trajectory.hpp"

#include <algorithm>
#include <list>
#include <map>
#include <vector>

#include "astronomy/epoch.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "numerics/fit_hermite_spline.hpp"

namespace principia {
namespace physics {
namespace internal_forkable {

using geometry::Instant;

template<typename Frame>
Instant const& ForkableTraits<DiscreteTrajectory<Frame>>::time(
    TimelineConstIterator const it) {
  return it->first;
}

template<typename Frame>
Instant const& DiscreteTrajectoryIterator<Frame>::time() const {
  return this->current()->first;
}

template<typename Frame>
DegreesOfFreedom<Frame> const&
DiscreteTrajectoryIterator<Frame>::degrees_of_freedom() const {
  return this->current()->second;
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::reference
DiscreteTrajectoryIterator<Frame>::operator*() const {
  auto const& it = this->current();
  return {it->first, it->second};
}

template<typename Frame>
std::optional<typename DiscreteTrajectoryIterator<Frame>::reference>
    DiscreteTrajectoryIterator<Frame>::operator->() const {
  auto const& it = this->current();
  return std::make_optional<reference>({it->first, it->second});
}

template<typename Frame>
not_null<DiscreteTrajectoryIterator<Frame>*>
DiscreteTrajectoryIterator<Frame>::that() {
  return this;
}

template<typename Frame>
not_null<DiscreteTrajectoryIterator<Frame> const*>
DiscreteTrajectoryIterator<Frame>::that() const {
  return this;
}

}  // namespace internal_forkable

namespace internal_discrete_trajectory {

using astronomy::InfiniteFuture;
using astronomy::InfinitePast;
using base::make_not_null_unique;
using numerics::FitHermiteSpline;

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkWithCopy(Instant const& time) {
  // May be at |timeline_end()| if |time| is the fork time of this object.
  auto timeline_it = timeline_.find(time);
  CHECK(timeline_it != timeline_end() ||
        (!this->is_root() && time == this->Fork().time()))
      << "NewForkWithCopy at nonexistent time " << time;

  auto const fork = this->NewFork(timeline_it);

  // Copy the tail of the trajectory in the child object.
  if (timeline_it != timeline_.end()) {
    fork->timeline_.insert(++timeline_it, timeline_.end());
  }
  return fork;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkWithoutCopy(Instant const& time) {
  // May be at |timeline_end()| if |time| is the fork time of this object.
  auto timeline_it = timeline_.find(time);
  CHECK(timeline_it != timeline_end() ||
        (!this->is_root() && time == this->Fork().time()))
      << "NewForkWithoutCopy at nonexistent time " << time;

  return this->NewFork(timeline_it);
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkAtLast() {
  auto end = timeline_.end();
  if (timeline_.empty()) {
    return this->NewFork(end);
  } else {
    return this->NewFork(--end);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::AttachFork(
    not_null<std::unique_ptr<DiscreteTrajectory<Frame>>> fork) {
  CHECK(fork->is_root());
  CHECK(!this->Empty());

  auto& fork_timeline = fork->timeline_;
  auto const this_last = --this->end();

  // Determine if |fork| already has a point matching the end of this
  // trajectory.
  bool must_prepend;
  if (fork_timeline.empty()) {
    must_prepend = true;
  } else {
    CHECK_LE(this_last->time, fork_timeline.begin()->first);
    auto const it = fork_timeline.find(this_last.time());
    if (it == fork_timeline.end()) {
      must_prepend = true;
    } else {
      CHECK(it == fork_timeline.begin())
          << it->first << " " << this_last.time();
      must_prepend = false;
    }
  }

  // If needed, prepend to |fork| a copy of the last point of this trajectory.
  // This ensures that |fork| and this trajectory start and end, respectively,
  // with points at the same time (but possibly distinct degrees of freedom).
  if (must_prepend) {
    fork_timeline.emplace_hint(fork_timeline.begin(),
                               this_last.time(),
                               this_last.degrees_of_freedom());
  }

  // Attach |fork| to this trajectory.
  this->AttachForkToCopiedBegin(std::move(fork));

  // Remove the first point of |fork| now that it is properly attached to its
  // parent: that point is either redundant (if it was prepended above) or wrong
  // (because we "trust" this trajectory more than |fork|).  The children that
  // might have been forked at the deleted point were relocated by
  // AttachForkToCopiedBegin.
  fork_timeline.erase(fork_timeline.begin());
}

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>>
DiscreteTrajectory<Frame>::DetachFork() {
  CHECK(!this->is_root());

  // Insert a new point in the timeline for the fork time.  It should go at the
  // beginning of the timeline.
  auto const fork_it = this->Fork();
  auto const begin_it = timeline_.emplace_hint(
      timeline_.begin(), fork_it.time(), fork_it.degrees_of_freedom());
  CHECK(begin_it == timeline_.begin());

  // Detach this trajectory and tell the caller that it owns the pieces.
  return this->DetachForkWithCopiedBegin();
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  CHECK(this->is_root() || time > this->Fork().time())
       << "Append at " << time << " which is before fork time "
       << this->Fork().time();

  if (!timeline_.empty() && timeline_.cbegin()->first == time) {
    LOG(WARNING) << "Append at existing time " << time
                 << ", time range = [" << this->begin().time() << ", "
                 << back().time << "]";
    return;
  }
  auto it = timeline_.emplace_hint(timeline_.end(),
                                   time,
                                   degrees_of_freedom);
  // Decrementing |end()| is much faster than incrementing |it|.  Don't ask.
  CHECK(--timeline_.end() == it)
      << "Append out of order at " << time << ", last time is "
      << (--timeline_.end())->first;
  if (downsampling_.has_value()) {
    if (timeline_.size() == 1) {
      downsampling_->SetStartOfDenseTimeline(timeline_.begin(), timeline_);
    } else {
      this->CheckNoForksBefore(back().time);
      downsampling_->increment_dense_intervals(timeline_);
      if (downsampling_->reached_max_dense_intervals()) {
        std::vector<TimelineConstIterator> dense_iterators;
        // This contains points, hence one more than intervals.
        dense_iterators.reserve(downsampling_->max_dense_intervals() + 1);
        for (TimelineConstIterator it =
                 downsampling_->start_of_dense_timeline();
             it != timeline_.end();
             ++it) {
          dense_iterators.push_back(it);
        }
        auto right_endpoints = FitHermiteSpline<Instant, Position<Frame>>(
            dense_iterators,
            [](auto&& it) -> auto&& { return it->first; },
            [](auto&& it) -> auto&& { return it->second.position(); },
            [](auto&& it) -> auto&& { return it->second.velocity(); },
            downsampling_->tolerance());
        if (right_endpoints.empty()) {
          right_endpoints.push_back(dense_iterators.end() - 1);
        }
        TimelineConstIterator left = downsampling_->start_of_dense_timeline();
        for (const auto& it_in_dense_iterators : right_endpoints) {
          TimelineConstIterator const right = *it_in_dense_iterators;
          timeline_.erase(++left, right);
          left = right;
        }
        downsampling_->SetStartOfDenseTimeline(left, timeline_);
      }
    }
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetAfter(Instant const& time) {
  this->DeleteAllForksAfter(time);

  // Get an iterator denoting the first entry with time > |time|.  Remove that
  // entry and all the entries that follow it.  This preserves any entry with
  // time == |time|.
  auto const first_removed_in_timeline = timeline_.upper_bound(time);
  Instant const* const first_removed_time =
      first_removed_in_timeline == timeline_.end()
          ? nullptr
          : &first_removed_in_timeline->first;
  if (downsampling_.has_value()) {
    if (first_removed_time != nullptr &&
        *first_removed_time <= downsampling_->first_dense_time()) {
      // The start of the dense timeline will be invalidated.
      if (first_removed_in_timeline == timeline_.begin()) {
        // The timeline will be empty after erasing.
        downsampling_->SetStartOfDenseTimeline(timeline_.end(), timeline_);
      } else {
        // Further points will be appended to the last remaining point, so this
        // is where the dense timeline will begin.
        auto last_kept_in_timeline = first_removed_in_timeline;
        --last_kept_in_timeline;
        downsampling_->SetStartOfDenseTimeline(last_kept_in_timeline,
                                               timeline_);
      }
    }
  }
  timeline_.erase(first_removed_in_timeline, timeline_.end());
  if (downsampling_.has_value()) {
    downsampling_->RecountDenseIntervals(timeline_);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetBefore(Instant const& time) {
  this->CheckNoForksBefore(time);

  // Get an iterator denoting the first entry with time >= |time|.  Remove all
  // the entries that precede it.  This preserves any entry with time == |time|.
  auto const first_kept_in_timeline = timeline_.lower_bound(time);
  if (downsampling_.has_value() &&
      (first_kept_in_timeline == timeline_.end() ||
       downsampling_->first_dense_time() < first_kept_in_timeline->first)) {
    // The start of the dense timeline will be invalidated.
    downsampling_->SetStartOfDenseTimeline(first_kept_in_timeline, timeline_);
  }
  timeline_.erase(timeline_.begin(), first_kept_in_timeline);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::SetDownsampling(
    std::int64_t const max_dense_intervals,
    Length const& tolerance) {
  CHECK(this->is_root());
  CHECK(!downsampling_.has_value());
  downsampling_.emplace(
      max_dense_intervals, tolerance, timeline_.begin(), timeline_);
}
template<typename Frame>
void DiscreteTrajectory<Frame>::ClearDownsampling() {
  downsampling_.reset();
}

template<typename Frame>
Instant DiscreteTrajectory<Frame>::t_min() const {
  return this->Empty() ? InfiniteFuture : this->begin().time();
}

template<typename Frame>
Instant DiscreteTrajectory<Frame>::t_max() const {
  return this->Empty() ? InfinitePast : back().time;
}

template<typename Frame>
Position<Frame> DiscreteTrajectory<Frame>::EvaluatePosition(
    Instant const& time) const {
  return GetInterpolation(time).Evaluate(time);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectory<Frame>::EvaluateVelocity(
    Instant const& time) const {;
  return GetInterpolation(time).EvaluateDerivative(time);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectory<Frame>::EvaluateDegreesOfFreedom(
    Instant const& time) const {
  auto const interpolation = GetInterpolation(time);
  return {interpolation.Evaluate(time), interpolation.EvaluateDerivative(time)};
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> const message,
    std::vector<DiscreteTrajectory<Frame>*> const& forks)
    const {
  CHECK(this->is_root());

  std::vector<DiscreteTrajectory<Frame>*> mutable_forks = forks;
  WriteSubTreeToMessage(message, mutable_forks);
  CHECK(std::all_of(mutable_forks.begin(),
                    mutable_forks.end(),
                    [](DiscreteTrajectory<Frame>* const fork) {
                      return fork == nullptr;
                    }));
}

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>>
DiscreteTrajectory<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<DiscreteTrajectory<Frame>**> const& forks) {
  auto trajectory = make_not_null_unique<DiscreteTrajectory>();
  CHECK(std::all_of(forks.begin(),
                    forks.end(),
                    [](DiscreteTrajectory<Frame>** const fork) {
                      return fork != nullptr && *fork == nullptr;
                    }));
  trajectory->FillSubTreeFromMessage(message, forks);
  return trajectory;
}

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
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_lower_bound(Instant const& time) const {
  return timeline_.lower_bound(time);
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::timeline_empty() const {
  return timeline_.empty();
}

template<typename Frame>
std::int64_t DiscreteTrajectory<Frame>::timeline_size() const {
  return timeline_.size();
}

template<typename Frame>
DiscreteTrajectory<Frame>::Downsampling::Downsampling(
    std::int64_t const max_dense_intervals,
    Length const tolerance,
    TimelineConstIterator const start_of_dense_timeline,
    Timeline const& timeline)
    : max_dense_intervals_(max_dense_intervals),
      tolerance_(tolerance),
      start_of_dense_timeline_(start_of_dense_timeline) {
  RecountDenseIntervals(timeline);
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::Downsampling::start_of_dense_timeline() const {
  return start_of_dense_timeline_;
}

template<typename Frame>
Instant const& DiscreteTrajectory<Frame>::Downsampling::first_dense_time()
    const {
  return start_of_dense_timeline_->first;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::SetStartOfDenseTimeline(
    TimelineConstIterator const value,
    Timeline const& timeline) {
  start_of_dense_timeline_ = value;
  RecountDenseIntervals(timeline);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::RecountDenseIntervals(
    Timeline const& timeline) {
  dense_intervals_ =
      std::distance(start_of_dense_timeline_, timeline.end()) - 1;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::increment_dense_intervals(
    Timeline const& timeline) {
  ++dense_intervals_;
  DCHECK_EQ(dense_intervals_,
            std::distance(start_of_dense_timeline_, timeline.end()) - 1);
}

template<typename Frame>
std::int64_t DiscreteTrajectory<Frame>::Downsampling::max_dense_intervals()
    const {
  return max_dense_intervals_;
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::Downsampling::reached_max_dense_intervals()
    const {
  return dense_intervals_ >= max_dense_intervals_;
}

template<typename Frame>
Length DiscreteTrajectory<Frame>::Downsampling::tolerance() const {
  return tolerance_;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::WriteToMessage(
    not_null<serialization::DiscreteTrajectory::Downsampling*> message,
    Timeline const& timeline) const {
  if (start_of_dense_timeline_ == timeline.end()) {
    message->clear_start_of_dense_timeline();
  } else {
    first_dense_time().WriteToMessage(
        message->mutable_start_of_dense_timeline());
  }
  message->set_max_dense_intervals(max_dense_intervals_);
  tolerance_.WriteToMessage(message->mutable_tolerance());
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::Downsampling
DiscreteTrajectory<Frame>::Downsampling::ReadFromMessage(
    serialization::DiscreteTrajectory::Downsampling const& message,
    Timeline const& timeline) {
  TimelineConstIterator start_of_dense_timeline;
  if (message.has_start_of_dense_timeline()) {
    start_of_dense_timeline = timeline.find(
        Instant::ReadFromMessage(message.start_of_dense_timeline()));
    CHECK(start_of_dense_timeline != timeline.end());
  } else {
    start_of_dense_timeline = timeline.end();
  }
  return Downsampling(message.max_dense_intervals(),
                      Length::ReadFromMessage(message.tolerance()),
                      start_of_dense_timeline,
                      timeline);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteSubTreeToMessage(
    not_null<serialization::DiscreteTrajectory*> const message,
    std::vector<DiscreteTrajectory<Frame>*>& forks) const {
  Forkable<DiscreteTrajectory, Iterator>::WriteSubTreeToMessage(message, forks);
  for (auto const& pair : timeline_) {
    Instant const& instant = pair.first;
    DegreesOfFreedom<Frame> const& degrees_of_freedom = pair.second;
    auto const instantaneous_degrees_of_freedom = message->add_timeline();
    instant.WriteToMessage(instantaneous_degrees_of_freedom->mutable_instant());
    degrees_of_freedom.WriteToMessage(
        instantaneous_degrees_of_freedom->mutable_degrees_of_freedom());
  }
  if (downsampling_.has_value()) {
    downsampling_->WriteToMessage(message->mutable_downsampling(), timeline_);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::FillSubTreeFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<DiscreteTrajectory<Frame>**> const& forks) {
  for (auto timeline_it = message.timeline().begin();
       timeline_it != message.timeline().end();
       ++timeline_it) {
    Append(Instant::ReadFromMessage(timeline_it->instant()),
           DegreesOfFreedom<Frame>::ReadFromMessage(
               timeline_it->degrees_of_freedom()));
  }
  if (message.has_downsampling()) {
    CHECK(this->is_root());
    downsampling_.emplace(
        Downsampling::ReadFromMessage(message.downsampling(), timeline_));
  }
  Forkable<DiscreteTrajectory, Iterator>::FillSubTreeFromMessage(message,
                                                                 forks);
}

template<typename Frame>
Hermite3<Instant, Position<Frame>> DiscreteTrajectory<Frame>::GetInterpolation(
    Instant const& time) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  // This is the upper bound of the interval upon which we will do the
  // interpolation.
  auto const upper = this->LowerBound(time);
  auto const lower = upper == this->begin() ? upper : --Iterator{upper};
  return Hermite3<Instant, Position<Frame>>{
      {lower.time(), upper.time()},
      {lower.degrees_of_freedom().position(),
       upper.degrees_of_freedom().position()},
      {lower.degrees_of_freedom().velocity(),
       upper.degrees_of_freedom().velocity()}};
}

}  // namespace internal_discrete_trajectory
}  // namespace physics
}  // namespace principia
