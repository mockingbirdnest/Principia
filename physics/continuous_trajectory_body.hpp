
#pragma once

#include <algorithm>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

#include "astronomy/epoch.hpp"
#include "glog/stl_logging.h"
#include "numerics/ulp_distance.hpp"
#include "physics/continuous_trajectory.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_continuous_trajectory {

using base::Error;
using base::make_not_null_unique;
using numerics::ULPDistance;
using quantities::DebugString;
using quantities::SIUnit;
using quantities::si::Metre;
using quantities::si::Second;

int const max_degree = 17;
int const min_degree = 3;
int const max_degree_age = 100;

// Only supports 8 divisions for now.
int const divisions = 8;

template<typename Frame>
ContinuousTrajectory<Frame>::ContinuousTrajectory(Time const& step,
                                                  Length const& tolerance)
    : step_(step),
      tolerance_(tolerance),
      adjusted_tolerance_(tolerance_),
      is_unstable_(false),
      degree_(min_degree),
      degree_age_(0) {
  CHECK_LT(0 * Metre, tolerance_);
}

template<typename Frame>
bool ContinuousTrajectory<Frame>::empty() const {
  return series_.empty();
}

template<typename Frame>
double ContinuousTrajectory<Frame>::average_degree() const {
  if (empty()) {
    return 0;
  } else {
    double total = 0;
    for (auto const& series : series_) {
      total += series.degree();
    }
    return total / series_.size();
  }
}

template<typename Frame>
Status ContinuousTrajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  // Consistency checks.
  if (first_time_) {
    Instant const t0;
    CHECK_GE(1,
             ULPDistance((last_points_.back().first + step_ - t0) /
                             SIUnit<Time>(),
                         (time - t0) / SIUnit<Time>()))
        << "Append at times that are not equally spaced, expected " << step_
        << ", found " << last_points_.back().first << " and " << time;
  } else {
    first_time_ = time;
  }

  Status status;
  if (last_points_.size() == divisions) {
    // These vectors are thread-local to avoid deallocation/reallocation each
    // time we go through this code path.
    thread_local std::vector<Displacement<Frame>> q(divisions + 1);
    thread_local std::vector<Velocity<Frame>> v(divisions + 1);
    q.clear();
    v.clear();

    for (auto const& pair : last_points_) {
      DegreesOfFreedom<Frame> const& degrees_of_freedom = pair.second;
      q.push_back(degrees_of_freedom.position() - Frame::origin);
      v.push_back(degrees_of_freedom.velocity());
    }
    q.push_back(degrees_of_freedom.position() - Frame::origin);
    v.push_back(degrees_of_freedom.velocity());

    status = ComputeBestNewhallApproximation(
        time, q, v, &ЧебышёвSeries<Displacement<Frame>>::NewhallApproximation);

    // Wipe-out the points that have just been incorporated in a series.
    last_points_.clear();
  }

  // Note that we only insert the new point in the map *after* computing the
  // approximation, because clearing the map is much more efficient than erasing
  // every element but one.
  last_points_.emplace_back(time, degrees_of_freedom);

  return status;
}

template<typename Frame>
void ContinuousTrajectory<Frame>::ForgetBefore(Instant const& time) {
  if (time < t_min()) {
    // TODO(phl): test for this case, it yielded a check failure in
    // |FindSeriesForInstant|.
    return;
  }
  series_.erase(series_.begin(), FindSeriesForInstant(time));

  // If there are no |series_| left, clear everything.  Otherwise, update the
  // first time.
  if (series_.empty()) {
    first_time_ = std::experimental::nullopt;
    last_points_.clear();
  } else {
    first_time_ = time;
  }
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_min() const {
  if (empty()) {
    return astronomy::InfiniteFuture;
  }
  return *first_time_;
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_max() const {
  if (empty()) {
    return astronomy::InfinitePast;
  }
  return series_.back().t_max();
}

template<typename Frame>
Position<Frame> ContinuousTrajectory<Frame>::EvaluatePosition(
    Instant const& time) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  auto const it = FindSeriesForInstant(time);
  CHECK(it != series_.end());
  return it->Evaluate(time) + Frame::origin;
}

template<typename Frame>
Velocity<Frame> ContinuousTrajectory<Frame>::EvaluateVelocity(
    Instant const& time) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  auto const it = FindSeriesForInstant(time);
  CHECK(it != series_.end());
  return it->EvaluateDerivative(time);
}

template<typename Frame>
DegreesOfFreedom<Frame> ContinuousTrajectory<Frame>::EvaluateDegreesOfFreedom(
    Instant const& time) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  auto const it = FindSeriesForInstant(time);
  CHECK(it != series_.end());
  return DegreesOfFreedom<Frame>(it->Evaluate(time) + Frame::origin,
                                 it->EvaluateDerivative(time));
}

template<typename Frame>
typename ContinuousTrajectory<Frame>::Checkpoint
ContinuousTrajectory<Frame>::GetCheckpoint() const {
  return {t_max(),
          adjusted_tolerance_,
          is_unstable_,
          degree_,
          degree_age_,
          last_points_};
}

template<typename Frame>
void ContinuousTrajectory<Frame>::WriteToMessage(
      not_null<serialization::ContinuousTrajectory*> const message) const {
  WriteToMessage(message, GetCheckpoint());
}

template<typename Frame>
void ContinuousTrajectory<Frame>::WriteToMessage(
      not_null<serialization::ContinuousTrajectory*> const message,
      Checkpoint const& checkpoint) const {
  step_.WriteToMessage(message->mutable_step());
  tolerance_.WriteToMessage(message->mutable_tolerance());
  checkpoint.adjusted_tolerance_.WriteToMessage(
      message->mutable_adjusted_tolerance());
  message->set_is_unstable(checkpoint.is_unstable_);
  message->set_degree(checkpoint.degree_);
  message->set_degree_age(checkpoint.degree_age_);
  for (auto const& s : series_) {
    if (s.t_max() <= checkpoint.t_max_) {
      s.WriteToMessage(message->add_series());
    }
    if (s.t_max() == checkpoint.t_max_) {
      break;
    }
    CHECK_LT(s.t_max(), checkpoint.t_max_);
  }
  if (first_time_) {
    first_time_->WriteToMessage(message->mutable_first_time());
  }
  for (auto const& pair : checkpoint.last_points_) {
    Instant const& instant = pair.first;
    DegreesOfFreedom<Frame> const& degrees_of_freedom = pair.second;
    not_null<
        serialization::ContinuousTrajectory::InstantaneousDegreesOfFreedom*>
        const instantaneous_degrees_of_freedom = message->add_last_point();
    instant.WriteToMessage(instantaneous_degrees_of_freedom->mutable_instant());
    degrees_of_freedom.WriteToMessage(
        instantaneous_degrees_of_freedom->mutable_degrees_of_freedom());
  }
}

template<typename Frame>
not_null<std::unique_ptr<ContinuousTrajectory<Frame>>>
ContinuousTrajectory<Frame>::ReadFromMessage(
      serialization::ContinuousTrajectory const& message) {
  not_null<std::unique_ptr<ContinuousTrajectory<Frame>>> continuous_trajectory =
      std::make_unique<ContinuousTrajectory<Frame>>(
          Time::ReadFromMessage(message.step()),
          Length::ReadFromMessage(message.tolerance()));
  continuous_trajectory->adjusted_tolerance_ =
      Length::ReadFromMessage(message.adjusted_tolerance());
  continuous_trajectory->is_unstable_ = message.is_unstable();
  continuous_trajectory->degree_ = message.degree();
  continuous_trajectory->degree_age_ = message.degree_age();
  for (auto const& s : message.series()) {
    continuous_trajectory->series_.push_back(
        ЧебышёвSeries<Displacement<Frame>>::ReadFromMessage(s));
  }
  if (message.has_first_time()) {
    continuous_trajectory->first_time_ =
        Instant::ReadFromMessage(message.first_time());
  }
  for (auto const& l : message.last_point()) {
    continuous_trajectory->last_points_.push_back(
        {Instant::ReadFromMessage(l.instant()),
         DegreesOfFreedom<Frame>::ReadFromMessage(l.degrees_of_freedom())});
  }
  return continuous_trajectory;
}

template<typename Frame>
bool ContinuousTrajectory<Frame>::Checkpoint::IsAfter(
    Instant const& time) const {
  return time < t_max_;
}

template<typename Frame>
ContinuousTrajectory<Frame>::Checkpoint::Checkpoint(
    Instant const& t_max,
    Length const& adjusted_tolerance,
    bool const is_unstable,
    int const degree,
    int const degree_age,
    std::vector<std::pair<Instant, DegreesOfFreedom<Frame>>> const& last_points)
    : t_max_(t_max),
      adjusted_tolerance_(adjusted_tolerance),
      is_unstable_(is_unstable),
      degree_(degree),
      degree_age_(degree_age),
      last_points_(last_points) {}

template<typename Frame>
ContinuousTrajectory<Frame>::ContinuousTrajectory() {}

template<typename Frame>
Status ContinuousTrajectory<Frame>::ComputeBestNewhallApproximation(
    Instant const& time,
    std::vector<Displacement<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v,
    ЧебышёвSeries<Displacement<Frame>> (*newhall_approximation)(
        int const degree,
        std::vector<Displacement<Frame>> const& q,
        std::vector<Velocity<Frame>> const& v,
        Instant const& t_min,
        Instant const& t_max)) {
  Length const previous_adjusted_tolerance = adjusted_tolerance_;

  // If the degree is too old, restart from the lowest degree.  This ensures
  // that we use the lowest possible degree at a small computational cost.
  if (degree_age_ >= max_degree_age) {
    VLOG(1) << "Lowering degree for " << this << " from " << degree_
            << " to " << min_degree << " because the approximation is too old";
    is_unstable_ = false;
    adjusted_tolerance_ = tolerance_;
    degree_ = min_degree;
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
    VLOG(1) << "Lowering degree for " << this << " from " << degree_
            << " to " << min_degree
            << " because error estimate " << error_estimate
            << " exceeds adjusted tolerance " << adjusted_tolerance_
            << " and computations are unstable";
    is_unstable_ = false;
    adjusted_tolerance_ = tolerance_;
    degree_ = min_degree - 1;
    degree_age_ = 0;
    previous_error_estimate = std::numeric_limits<double>::max() * Metre;
    error_estimate = 0.5 * previous_error_estimate;
  }

  // Increase the degree if the approximation is not accurate enough.  Stop
  // when we reach the maximum degree or when the error estimate is not
  // decreasing.
  while (error_estimate > adjusted_tolerance_ &&
         error_estimate < previous_error_estimate &&
         degree_ < max_degree) {
    ++degree_;
    VLOG(1) << "Increasing degree for " << this << " to " <<degree_
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
    if (degree_ > min_degree) {
      --degree_;
    }
    VLOG(1) << "Reverting to degree " << degree_ << " for " << this
            << " because error estimate increased (" << error_estimate
            << " vs. " << previous_error_estimate << ")";
    is_unstable_ = true;
    error_estimate = previous_error_estimate;
    adjusted_tolerance_ = std::max(adjusted_tolerance_, error_estimate);
  } else {
    VLOG(1) << "Using degree " << degree_ << " for " << this
            << " with error estimate " << error_estimate;
  }

  ++degree_age_;

  // Check that the tolerance did not explode.
  if (adjusted_tolerance_ < 1e6 * previous_adjusted_tolerance) {
    return Status::OK;
  } else {
    std::stringstream message;
    message << "Error trying to fit a smooth polynomial to the trajectory. "
            << "The approximation error jumped from "
            << previous_adjusted_tolerance << " to " << adjusted_tolerance_
            << " at time " << time << ". The last position is " << q.back()
            << " and the last velocity is " << v.back()
            << ". An apocalypse occurred and two celestials probably "
            << "collided because your solar system is unstable.";
    return Status(Error::INVALID_ARGUMENT, message.str());
  }
}

template<typename Frame>
typename std::vector<ЧебышёвSeries<Displacement<Frame>>>::const_iterator
ContinuousTrajectory<Frame>::FindSeriesForInstant(Instant const& time) const {
  // Need to use |lower_bound|, not |upper_bound|, because it allows
  // heterogeneous arguments.  This returns the first series |s| such that
  // |time <= s.t_max()|.
  auto const it = std::lower_bound(
                      series_.begin(), series_.end(), time,
                      [](ЧебышёвSeries<Displacement<Frame>> const& left,
                         Instant const& right) {
                        return left.t_max() < right;
                      });
  return it;
}

}  // namespace internal_continuous_trajectory
}  // namespace physics
}  // namespace principia
