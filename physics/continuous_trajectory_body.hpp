#pragma once

#include <algorithm>

#include "physics/continuous_trajectory.hpp"

namespace principia {
namespace physics {

namespace {

int const kMaxDegree = 17;
int const kMinDegree = 3;

// Only supports 8 divisions for now.
int const kDivisions = 8;

}  // namespace

template<typename Frame>
ContinuousTrajectory<Frame>::ContinuousTrajectory(Time const& step,
                                                  Length const& low_tolerance,
                                                  Length const& high_tolerance)
    : step_(step),
      low_tolerance_(low_tolerance),
      high_tolerance_(high_tolerance),
      degree_((kMinDegree + kMaxDegree) / 2) {
  CHECK_LT(low_tolerance_, high_tolerance_);
}

template<typename Frame>
bool ContinuousTrajectory<Frame>::empty() const {
  return series_.empty();
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_min() const {
  CHECK(!empty()) << "Empty trajectory";
  return series_.front().t_min();
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_max() const {
  CHECK(!empty()) << "Empty trajectory";
  return series_.back().t_max();
}

template<typename Frame>
void ContinuousTrajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  // Consistency checks.
  if (empty()) {
    first_time_ = std::make_unique<Instant>(time);
  } else {
    CHECK_EQ(last_points_.back().first + step_, time)
        << "Append at times that are not equally spaced";
  }

  if (last_points_.size() == kDivisions) {
    // These vectors are static to avoid deallocation/reallocation each time we
    // go through this code path.
    static std::vector<Displacement<Frame>> q(kDivisions + 1);
    static std::vector<Velocity<Frame>> v(kDivisions + 1);
    q.clear();
    v.clear();

    for (auto const& pair : last_points_) {
      DegreesOfFreedom<Frame> const& degrees_of_freedom = pair.second;
      q.push_back(degrees_of_freedom.position() - Frame::origin);
      v.push_back(degrees_of_freedom.velocity());
    }
    q.push_back(degrees_of_freedom.position() - Frame::origin);
    v.push_back(degrees_of_freedom.velocity());
    LOG(ERROR)<<q.size()<<" "<<v.size()<<" "<<degree_;

    // Compute the approximation with the current degree.
    series_.push_back(
        ЧебышёвSeries<Displacement<Frame>>::NewhallApproximation(
            degree_, q, v, last_points_.cbegin()->first, time));

    Length error_estimate = series_.back().last_coefficient().Norm();

    // Increase the degree if the approximation is not accurate enough.
    while (error_estimate > high_tolerance_ && degree_ < kMaxDegree) {
      ++degree_;
      series_.back() =
          ЧебышёвSeries<Displacement<Frame>>::NewhallApproximation(
              degree_, q, v, last_points_.cbegin()->first, time);

      error_estimate = series_.back().last_coefficient().Norm();
    }

    // Try to decrease the degree if the approximation is too accurate, but make
    // sure that we don't go above |high_tolerance_|.
    while (error_estimate < low_tolerance_ && degree_ > kMinDegree) {
      --degree_;
      auto tentative_series =
          ЧебышёвSeries<Displacement<Frame>>::NewhallApproximation(
              degree_, q, v, last_points_.cbegin()->first, time);

      error_estimate = series_.back().last_coefficient().Norm();
      if (error_estimate > high_tolerance_) {
        ++degree_;
      } else {
        series_.back() = std::move(tentative_series);
      }
    }

    // Wipe-out the vector.
    last_points_.clear();
  }

  // Note that we only insert the new point in the map *after* computing the
  // approximation, because clearing the map is much more efficient than erasing
  // every element but one.
  last_points_.emplace_back(time, degrees_of_freedom);
}

template<typename Frame>
void ContinuousTrajectory<Frame>::ForgetBefore(Instant const& time) {
  series_.erase(series_.begin(), FindSeriesForInstant(time));

  // If there are no |series_| left, clear everything.  Otherwise, update the
  // first time.
  if (series_.empty()) {
    first_time_.reset();
    last_points_.clear();
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
typename std::vector<ЧебышёвSeries<Displacement<Frame>>>::const_iterator
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
