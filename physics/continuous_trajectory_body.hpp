#pragma once

#include <algorithm>

#include "physics/continuous_trajectory.hpp"

namespace principia {
namespace physics {

template<typename Frame>
ContinuousTrajectory<Frame>::ContinuousTrajectory(int const degree)
    : degree_(degree) {}

template<typename Frame>
bool ContinuousTrajectory<Frame>::empty() const {
  return series_.empty();
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::first_time() const {
  CHECK(!empty()) << "Empty trajectory";
  return series_.front().t_min();
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::last_time() const {
  CHECK(!empty()) << "Empty trajectory";
  return series_.back().t_max();
}

template<typename Frame>
void ContinuousTrajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  // Only supports 8 divisions for now.
  int const kDivisions = 8;

  // Consistency checks.
  if (empty()) {
    first_time_ = std::make_unique<Instant>(time);
  } else {
    if (interval_ == nullptr) {
      interval_ = std::make_unique<Time>(time - *first_time_);
      CHECK_GT(Time{}, *interval_)
          << "Append before the end of the trajectory";
    } else {
      CHECK_EQ(*interval_, time - last_points_.crbegin()->first)
          << "Append at times that are not equally spaced";
    }
  }

  if (last_points_.size() == kDivisions) {
    // These vectors are static to avoid deallocation/reallocation each time we
    // go through this code path.
    static std::vector<Length> q(kDivisions + 1);
    static std::vector<Speed> v(kDivisions + 1);
    q.clear();
    v.clear();

    for (auto const& pair : last_points_) {
      DegreesOfFreedom<Frame> const& degrees_of_freedom = pair.second;
      q.push_back(degrees_of_freedom.position - Frame::origin);
      v.push_back(degrees_of_freedom.velocity);
    }
    q.push_back(degrees_of_freedom.position - Frame::origin);
    v.push_back(degrees_of_freedom.velocity);

    // Compute the approximation.
    series_.push_back(ЧебышёвSeries::NewhallApproximation(
        degree_, q, v, last_points_.cbegin()->first, time);

    // Wipe-out the map.
    last_points_.clear();
  }

  // Note that we only insert the new point in the map *after* computing the
  // approximation, because clearing the map is much more efficient than erasing
  // every element but one.
  last_points[time] = degrees_of_freedom;
}

template<typename Frame>
void ContinuousTrajectory<Frame>::ForgetBefore(Instant const& time) {
  series_.erase(series_.begin, FindSeriesForInstant(time));

  // If there are no |series_| left, clear everything.  Otherwise, update the
  // first time.
  if (series_.empty()) {
    interval_.reset();
    first_time_.reset();
    last_point_.reset();
  } else {
    *first_time_ = time;
  }
}

template<typename Frame>
Position<Frame> ContinuousTrajectory<Frame>::EvaluatePosition(
    Instant const& time,
    Hint* const hint) const {
  if (MayUseHint(time, hint)) {
    return series_[index].Evaluate(time);
  } else {
    auto const it = FindSeriesForInstant(time);
    if (hint != nullptr) {
      *hint = it - series_.cbegin();
    }
    return it->Evaluate(time);
  }
}

template<typename Frame>
Velocity<Frame> ContinuousTrajectory<Frame>::EvaluateVelocity(
    Instant const& time,
    Hint* const hint) const {
  if (MayUseHint(time, hint)) {
    return series_[index].EvaluateDerivative(time);
  } else {
    auto const it = FindSeriesForInstant(time);
    if (hint != nullptr) {
      *hint = it - series_.cbegin();
    }
    return it->EvaluateDerivative(time);
  }
}

template<typename Frame>
DegreesOfFreedom<Frame> ContinuousTrajectory<Frame>::EvaluateDegreesOfFreedom(
    Instant const& time,
    Hint* const hint) const {
  if (MayUseHint(time, hint)) {
    ЧебышёвSeries const& series = series_[index];
    return DegreesOfFreedom<Frame>(series.Evaluate(time),
                                   series.EvaluateDerivative(time));
  } else {
    auto const it = FindSeriesForInstant(time);
    if (hint != nullptr) {
      *hint = it - series_.cbegin();
    }
    return DegreesOfFreedom<Frame>(it->Evaluate(time),
                                   it->EvaluateDerivative(time));
  }
}

template<typename Frame>
ContinuousTrajectory<Frame>::Hint::Hint()
    : index_(std::numeric_limits<int>::max()) {}

template<typename Frame>
std::vector<ЧебышёвSeries>::const_iterator
ContinuousTrajectory<Frame>::FindSeriesForInstant(Instant const& time) const {
  return std::upper_bound(series_.begin(), series_.end(), time,
                          [](ЧебышёвSeries const& left, Instant const& right) {
                            return left.t_max() < right;
                          });
}

template<typename Frame>
bool ContinuousTrajectory<Frame>::MayUseHint(Instant const& time,
                                             Hint* const hint) const {
  if (hint != nullptr) {
    // A shorthand for the index held by the |hint|.
    int& index = hint->index_;
    if (index < series_.size() && series_[index].t_min() <= time) {
      if (time <= series_[index].t_max()) {
        // Use this interval.
        return true;
      } else if (index < series_.size() - 1 &&
                 time <= series_[index + 1].t_max()) {
        // Move to the next interval.
        ++index;
        return true;
      }
    }
  }
  return false;
}

}  // namespace physics
}  // namespace principia
