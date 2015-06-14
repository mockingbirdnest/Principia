#pragma once

#include <algorithm>
#include <limits>
#include <vector>

#include "physics/continuous_trajectory.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using testing_utilities::ULPDistance;

namespace physics {

namespace {

int const kMaxDegree = 17;
int const kMinDegree = 3;
int const kMaxDegreeAge = 100;

// Only supports 8 divisions for now.
int const kDivisions = 8;

}  // namespace

template<typename Frame>
ContinuousTrajectory<Frame>::ContinuousTrajectory(Time const& step,
                                                  Length const& low_tolerance,
                                                  Length const& high_tolerance)
    : step_(step),
      tolerance_(high_tolerance),
      adjusted_tolerance_(tolerance_),
      is_unstable_(false),
      degree_(kMinDegree),
      degree_age_(0) {
  CHECK_LT(0 * Metre, tolerance_);
}

template<typename Frame>
bool ContinuousTrajectory<Frame>::empty() const {
  return series_.empty();
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_min() const {
  CHECK(!empty()) << "Empty trajectory";
  return *first_time_;
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
  if (first_time_ == nullptr) {
    first_time_ = std::make_unique<Instant>(time);
  } else {
    Instant const t0;
    CHECK_GE(1,
             ULPDistance((last_points_.back().first + step_ - t0) /
                             SIUnit<Time>(),
                         (time - t0) / SIUnit<Time>()))
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

    ComputeBestNewhallApproximation(
        time, q, v, &ЧебышёвSeries<Displacement<Frame>>::NewhallApproximation);

    // Wipe-out the points that have just been incorporated in a series.
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
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  if (MayUseHint(time, hint)) {
    return series_[hint->index_].Evaluate(time) + Frame::origin;
  } else {
    auto const it = FindSeriesForInstant(time);
    if (hint != nullptr) {
      hint->index_ = it - series_.cbegin();
    }
    return it->Evaluate(time) + Frame::origin;
  }
}

template<typename Frame>
Velocity<Frame> ContinuousTrajectory<Frame>::EvaluateVelocity(
    Instant const& time,
    Hint* const hint) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  if (MayUseHint(time, hint)) {
    return series_[hint->index_].EvaluateDerivative(time);
  } else {
    auto const it = FindSeriesForInstant(time);
    if (hint != nullptr) {
      hint->index_ = it - series_.cbegin();
    }
    return it->EvaluateDerivative(time);
  }
}

template<typename Frame>
DegreesOfFreedom<Frame> ContinuousTrajectory<Frame>::EvaluateDegreesOfFreedom(
    Instant const& time,
    Hint* const hint) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  if (MayUseHint(time, hint)) {
    ЧебышёвSeries<Displacement<Frame>> const& series = series_[hint->index_];
    return DegreesOfFreedom<Frame>(series.Evaluate(time) + Frame::origin,
                                   series.EvaluateDerivative(time));
  } else {
    auto const it = FindSeriesForInstant(time);
    if (hint != nullptr) {
      hint->index_ = it - series_.cbegin();
    }
    return DegreesOfFreedom<Frame>(it->Evaluate(time) + Frame::origin,
                                   it->EvaluateDerivative(time));
  }
}

template<typename Frame>
ContinuousTrajectory<Frame>::Hint::Hint()
    : index_(std::numeric_limits<int>::max()) {}

  
template<typename Frame>
void ContinuousTrajectory<Frame>::ComputeBestNewhallApproximation(
    Instant const& time,
    std::vector<Displacement<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v,
    ЧебышёвSeries<Displacement<Frame>> (*newhall_approximation)(
        int const degree,
        std::vector<Displacement<Frame>> const& q,
        std::vector<Velocity<Frame>> const& v,
        Instant const& t_min,
        Instant const& t_max)) {
  // If the degree is too old, restart from the lowest degree.  This ensures
  // that we use the lowest possible degree at a small computational cost.
  LOG(ERROR)<<degree_age_;
  if (degree_age_ > kMaxDegreeAge) {
    LOG(ERROR) << "Lowering degree from " << degree_ << " to " << kMinDegree
            << " because the approximation is too old";
    is_unstable_ = false;
    adjusted_tolerance_ = tolerance_;
    degree_ = kMinDegree;
    degree_age_ = 0;
  }

  // Compute the approximation with the current degree.
  series_.push_back(
      newhall_approximation(degree_, q, v, last_points_.cbegin()->first, time));

  // Estimate the error.  For initializing |previous_error_estimate|, any value
  // greater than |error_estimate| will do.
  Length error_estimate = series_.back().last_coefficient().Norm();
  Length previous_error_estimate = error_estimate + error_estimate;

  // If we are in the zone of numerical instabilities and we exceeded the
  // tolerance, restart from the lowest degree.
  if (is_unstable_ && error_estimate > adjusted_tolerance_) {
    LOG(ERROR) << "Lowering degree from " << degree_ << " to " << kMinDegree
            << " because error estimate " << error_estimate
            << " exceeds adjusted tolerance " << adjusted_tolerance_
            << " and computations are unstable";
    is_unstable_ = false;
    adjusted_tolerance_ = tolerance_;
    degree_ = kMinDegree - 1;
    degree_age_ = 0;
    previous_error_estimate = std::numeric_limits<double>::max() * Metre;
    error_estimate = 0.5 * previous_error_estimate;
  }

  // Increase the degree if the approximation is not accurate enough.  Stop
  // when we reach the maximum degree or when the error estimate is not
  // decreasing.
  while (error_estimate > adjusted_tolerance_ &&
         error_estimate < previous_error_estimate &&
         degree_ < kMaxDegree) {
    ++degree_;
    LOG(ERROR) << "Increasing degree for " << this << " to " <<degree_
            << " because error estimate was " << error_estimate;
    series_.back() =
        newhall_approximation(
            degree_, q, v, last_points_.cbegin()->first, time);
    previous_error_estimate = error_estimate;
    error_estimate = series_.back().last_coefficient().Norm();
  }

  // If we have entered the zone of numerical instability, go back to the
  // point where the error was decreasing and nudge the tolerance since we
  // won't be able to reliably do better than that.
  if (error_estimate >= previous_error_estimate) {
    if (degree_ > kMinDegree) {
    --degree_;
    }
    LOG(ERROR) << "Reverting to degree " << degree_ << " for " << this
            << " because error estimate increased (" << error_estimate
            << " vs. " << previous_error_estimate << ")";
    is_unstable_ = true;
    error_estimate = previous_error_estimate;
    adjusted_tolerance_ = std::max(adjusted_tolerance_, error_estimate);
  } else {
    LOG(ERROR) << "Using degree " << degree_ << " for " << this
            << " with error estimate " << error_estimate;
  }

  ++degree_age_;
}

template<typename Frame>
typename std::vector<ЧебышёвSeries<Displacement<Frame>>>::const_iterator
ContinuousTrajectory<Frame>::FindSeriesForInstant(Instant const& time) const {
  // Need to use |lower_bound|, not |upper_bound|, because it allows
  // heterogeneous arguments.
  auto const it = std::lower_bound(series_.begin(), series_.end(), time,
                      [](ЧебышёвSeries<Displacement<Frame>> const& left,
                         Instant const& right) {
                        return left.t_max() < right;
                      });
  CHECK(it != series_.end());
  return it;
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
