
#pragma once

#include "physics/continuous_trajectory.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <sstream>
#include <utility>
#include <vector>

#include "astronomy/epoch.hpp"
#include "glog/stl_logging.h"
#include "numerics/newhall.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/ulp_distance.hpp"
#include "numerics/чебышёв_series.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_continuous_trajectory {

using base::Error;
using base::make_not_null_unique;
using numerics::EstrinEvaluator;
using numerics::ULPDistance;
using numerics::ЧебышёвSeries;
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
  return polynomials_.empty();
}

template<typename Frame>
double ContinuousTrajectory<Frame>::average_degree() const {
  if (polynomials_.empty()) {
    return 0;
  } else {
    double total = 0;
    for (auto const& pair : polynomials_) {
      total += pair.polynomial->degree();
    }
    return total / polynomials_.size();
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

    status = ComputeBestNewhallApproximation(time, q, v);

    // Wipe-out the points that have just been incorporated in a polynomial.
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
    // |FindPolynomialForInstant|.
    return;
  }
  polynomials_.erase(polynomials_.begin(), FindPolynomialForInstant(time));

  // If there are no |polynomials_| left, clear everything.  Otherwise, update
  // the first time.
  if (polynomials_.empty()) {
    first_time_ = std::nullopt;
    last_points_.clear();
    last_accessed_polynomial_ = 0;
  } else {
    first_time_ = time;
    last_accessed_polynomial_ = polynomials_.size() - 1;
  }
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_min() const {
  if (polynomials_.empty()) {
    return astronomy::InfiniteFuture;
  }
  return *first_time_;
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_max() const {
  if (polynomials_.empty()) {
    return astronomy::InfinitePast;
  }
  return polynomials_.crbegin()->t_max;
}

template<typename Frame>
Position<Frame> ContinuousTrajectory<Frame>::EvaluatePosition(
    Instant const& time) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  auto const it = FindPolynomialForInstant(time);
  CHECK(it != polynomials_.end());
  auto const& polynomial = it->polynomial;
  return polynomial->Evaluate(time) + Frame::origin;
}

template<typename Frame>
Velocity<Frame> ContinuousTrajectory<Frame>::EvaluateVelocity(
    Instant const& time) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  auto const it = FindPolynomialForInstant(time);
  CHECK(it != polynomials_.end());
  auto const& polynomial = it->polynomial;
  return polynomial->EvaluateDerivative(time);
}

template<typename Frame>
DegreesOfFreedom<Frame> ContinuousTrajectory<Frame>::EvaluateDegreesOfFreedom(
    Instant const& time) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  auto const it = FindPolynomialForInstant(time);
  CHECK(it != polynomials_.end());
  auto const& polynomial = it->polynomial;
  return DegreesOfFreedom<Frame>(polynomial->Evaluate(time) + Frame::origin,
                                 polynomial->EvaluateDerivative(time));
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
  for (auto const& pair : polynomials_) {
    Instant const& t_max = pair.t_max;
    auto const& polynomial = pair.polynomial;
    if (t_max <= checkpoint.t_max_) {
      auto* const pair = message->add_instant_polynomial_pair();
      t_max.WriteToMessage(pair->mutable_t_max());
      polynomial->WriteToMessage(pair->mutable_polynomial());
    }
    if (t_max == checkpoint.t_max_) {
      break;
    }
    CHECK_LT(t_max, checkpoint.t_max_);
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
  bool const is_pre_cohen = message.series_size() > 0;
  not_null<std::unique_ptr<ContinuousTrajectory<Frame>>> continuous_trajectory =
      std::make_unique<ContinuousTrajectory<Frame>>(
          Time::ReadFromMessage(message.step()),
          Length::ReadFromMessage(message.tolerance()));
  continuous_trajectory->adjusted_tolerance_ =
      Length::ReadFromMessage(message.adjusted_tolerance());
  continuous_trajectory->is_unstable_ = message.is_unstable();
  continuous_trajectory->degree_ = message.degree();
  continuous_trajectory->degree_age_ = message.degree_age();
  if (is_pre_cohen) {
    for (auto const& s : message.series()) {
      // Read the series, evaluate it and use the resulting values to build a
      // polynomial in the monomial basis.
      auto const series =
          ЧебышёвSeries<Displacement<Frame>>::ReadFromMessage(s);
      Time const step = (series.t_max() - series.t_min()) / divisions;
      Instant t = series.t_min();
      std::vector<Displacement<Frame>> q;
      std::vector<Velocity<Frame>> v;
      for (int i = 0; i <= divisions; t += step, ++i) {
        q.push_back(series.Evaluate(t));
        v.push_back(series.EvaluateDerivative(t));
      }
      Displacement<Frame> error_estimate;  // Should we do something with this?
      continuous_trajectory->polynomials_.emplace_back(
          series.t_max(),
          continuous_trajectory->NewhallApproximationInMonomialBasis(
              series.degree(),
              q, v,
              series.t_min(), series.t_max(),
              error_estimate));
    }
  } else {
    for (auto const& pair : message.instant_polynomial_pair()) {
      continuous_trajectory->polynomials_.emplace_back(
          Instant::ReadFromMessage(pair.t_max()),
          Polynomial<Displacement<Frame>, Instant>::template ReadFromMessage<
              EstrinEvaluator>(pair.polynomial()));
    }
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
ContinuousTrajectory<Frame>::InstantPolynomialPair::InstantPolynomialPair(
    Instant const t_max,
    not_null<std::unique_ptr<Polynomial<Displacement<Frame>, Instant>>>
        polynomial)
    : t_max(t_max),
      polynomial(std::move(polynomial)) {}

template<typename Frame>
not_null<std::unique_ptr<Polynomial<Displacement<Frame>, Instant>>>
ContinuousTrajectory<Frame>::NewhallApproximationInMonomialBasis(
    int degree,
    std::vector<Displacement<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v,
    Instant const& t_min,
    Instant const& t_max,
    Displacement<Frame>& error_estimate) const {
  return numerics::NewhallApproximationInMonomialBasis<
            Displacement<Frame>, EstrinEvaluator>(degree,
                                                  q, v,
                                                  t_min, t_max,
                                                  error_estimate);
}

template<typename Frame>
Status ContinuousTrajectory<Frame>::ComputeBestNewhallApproximation(
    Instant const& time,
    std::vector<Displacement<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v) {
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
  Displacement<Frame> displacement_error_estimate;
  polynomials_.emplace_back(time,
                            NewhallApproximationInMonomialBasis(
                                degree_,
                                q, v,
                                last_points_.cbegin()->first, time,
                                displacement_error_estimate));

  // Estimate the error.  For initializing |previous_error_estimate|, any value
  // greater than |error_estimate| will do.
  Length error_estimate = displacement_error_estimate.Norm();
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
    polynomials_.back().polynomial = NewhallApproximationInMonomialBasis(
                                         degree_,
                                         q, v,
                                         last_points_.cbegin()->first, time,
                                         displacement_error_estimate);
    previous_error_estimate = error_estimate;
    error_estimate = displacement_error_estimate.Norm();
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
typename ContinuousTrajectory<Frame>::InstantPolynomialPairs::const_iterator
ContinuousTrajectory<Frame>::FindPolynomialForInstant(
    Instant const& time) const {
  // This returns the first polynomial |p| such that |time <= p.t_max|.
  {
    auto const begin = polynomials_.begin();
    auto const it = begin + last_accessed_polynomial_;
    if (it != polynomials_.end() && time <= it->t_max &&
        (it == begin || std::prev(it)->t_max < time)) {
      return it;
    }
  }
  {
    auto const it =
        std::lower_bound(polynomials_.begin(),
                         polynomials_.end(),
                         time,
                         [](InstantPolynomialPair const& left,
                            Instant const& right) {
                           return left.t_max < right;
                         });
    last_accessed_polynomial_ = it - polynomials_.begin();
    return it;
  }
}

}  // namespace internal_continuous_trajectory
}  // namespace physics
}  // namespace principia
