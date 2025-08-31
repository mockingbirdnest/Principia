#pragma once

#include "physics/continuous_trajectory.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <utility>
#include <vector>

#include "geometry/interval.hpp"
#include "glog/stl_logging.h"
#include "numerics/fma.hpp"
#include "numerics/newhall.hpp"
#include "numerics/poisson_series.hpp"
#include "numerics/polynomial_in_чебышёв_basis.hpp"
#include "numerics/ulp_distance.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _continuous_trajectory {
namespace internal {

using namespace principia::geometry::_interval;
using namespace principia::numerics::_fma;
using namespace principia::numerics::_newhall;
using namespace principia::numerics::_poisson_series;
using namespace principia::numerics::_polynomial_in_чебышёв_basis;
using namespace principia::numerics::_ulp_distance;
using namespace principia::quantities::_si;

constexpr int max_degree = 17;
constexpr int min_degree = 3;
int const max_degree_age = 100;

// Only supports 8 divisions for now.
int const divisions = 8;

template<typename Frame>
ContinuousTrajectory<Frame>::ContinuousTrajectory(Time const& step,
                                                  Length const& tolerance)
    : step_(step),
      tolerance_(tolerance),
      checkpointer_(
          make_not_null_unique<
              Checkpointer<serialization::ContinuousTrajectory>>(
                  MakeCheckpointerWriter(),
                  MakeCheckpointerReader())),
      adjusted_tolerance_(tolerance_),
      is_unstable_(false),
      degree_(min_degree),
      degree_age_(0),
      polynomial_evaluator_policy_(CanUseHardwareFMA
                                       ? Policy::AlwaysEstrin()
                                       : Policy::AlwaysEstrinWithoutFMA()) {
  CHECK_LT(0 * Metre, tolerance_);
}

template<typename Frame>
bool ContinuousTrajectory<Frame>::empty() const {
  absl::ReaderMutexLock l(&lock_);
  return polynomials_.empty();
}

template<typename Frame>
double ContinuousTrajectory<Frame>::average_degree() const {
  absl::ReaderMutexLock l(&lock_);
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
absl::Status ContinuousTrajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  absl::MutexLock l(&lock_);

  // Consistency checks.
  if (first_time_) {
    Instant const t0;
    CHECK_GE(1,
             ULPDistance((last_points_.back().first + step_ - t0) /
                             si::Unit<Time>,
                         (time - t0) / si::Unit<Time>))
        << "Append at times that are not equally spaced, expected " << step_
        << ", found " << last_points_.back().first << " and " << time;
  } else {
    first_time_ = time;
  }

  absl::Status status;
  CHECK_LE(last_points_.size(), divisions);
  if (last_points_.size() == divisions) {
    // These vectors are thread-local to avoid deallocation/reallocation each
    // time we go through this code path.
    thread_local std::vector<Position<Frame>> q(divisions + 1);
    thread_local std::vector<Velocity<Frame>> v(divisions + 1);
    q.clear();
    v.clear();

    for (auto const& [_, degrees_of_freedom] : last_points_) {
      q.push_back(degrees_of_freedom.position());
      v.push_back(degrees_of_freedom.velocity());
    }
    q.push_back(degrees_of_freedom.position());
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
void ContinuousTrajectory<Frame>::Prepend(ContinuousTrajectory&& prefix) {
  absl::MutexLock l1(&lock_);
  absl::MutexLock l2(&prefix.lock_);

  CHECK_EQ(step_, prefix.step_);
  CHECK_EQ(tolerance_, prefix.tolerance_);

  if (prefix.polynomials_.empty()) {
    // Nothing to do.
  } else if (polynomials_.empty()) {
    // All the data comes from `prefix`.  This must set all the fields of
    // this object that are not set at construction.
    adjusted_tolerance_ = prefix.adjusted_tolerance_;
    is_unstable_ = prefix.is_unstable_;
    degree_ = prefix.degree_;
    degree_age_ = prefix.degree_age_;
    polynomials_ = std::move(prefix.polynomials_);
    last_accessed_polynomial_ = prefix.last_accessed_polynomial_;
    first_time_ = prefix.first_time_;
    last_points_ = prefix.last_points_;
  } else {
    // The polynomials must be aligned, because the time computations only use
    // basic arithmetic and are platform-independent.  The space computations,
    // on the other may depend on characteristics of the hardware and/or math
    // library, so we cannot check that the trajectories are "continuous" at the
    // junction.
    CHECK_EQ(*first_time_, prefix.polynomials_.back().t_max);
    // This operation is in O(prefix.size()).
    std::move(polynomials_.begin(),
              polynomials_.end(),
              std::back_inserter(prefix.polynomials_));
    polynomials_.swap(prefix.polynomials_);
    first_time_ = prefix.first_time_;
    // Note that any `last_points_` in `prefix` are irrelevant because they
    // correspond to a time interval covered by the first polynomial of this
    // object.
  }
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_min() const {
  absl::ReaderMutexLock l(&lock_);
  return t_min_locked();
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_max() const {
  absl::ReaderMutexLock l(&lock_);
  return t_max_locked();
}

template<typename Frame>
Position<Frame> ContinuousTrajectory<Frame>::EvaluatePosition(
    Instant const& time) const {
  absl::ReaderMutexLock l(&lock_);
  return EvaluatePositionLocked(time);
}

template<typename Frame>
Velocity<Frame> ContinuousTrajectory<Frame>::EvaluateVelocity(
    Instant const& time) const {
  absl::ReaderMutexLock l(&lock_);
  return EvaluateVelocityLocked(time);
}

template<typename Frame>
DegreesOfFreedom<Frame> ContinuousTrajectory<Frame>::EvaluateDegreesOfFreedom(
    Instant const& time) const {
  absl::ReaderMutexLock l(&lock_);
  return EvaluateDegreesOfFreedomLocked(time);
}

#if PRINCIPIA_CONTINUOUS_TRAJECTORY_SUPPORTS_PIECEWISE_POISSON_SERIES

template<typename Frame>
int ContinuousTrajectory<Frame>::PiecewisePoissonSeriesDegree(
    Instant const& t_min,
    Instant const& t_max) const {
  absl::ReaderMutexLock l(&lock_);
  CHECK_LE(t_min_locked(), t_min);
  CHECK_GE(t_max_locked(), t_max);
  auto const it_min = FindPolynomialForInstantLocked(t_min);
  auto const it_max = FindPolynomialForInstantLocked(t_max);
  int degree = min_degree;
  for (auto it = it_min;; ++it) {
    degree = std::max(degree, it->polynomial->degree());
    if (it == it_max) {
      break;
    }
  }
  return degree;
}


// Casts the polynomial to one of degree d1, and then increases the degree to
// d2.
#define PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(polynomial, d1, d2)     \
  case d1: {                                                                   \
    if constexpr (d1 <= d2) {                                                  \
      return PolynomialInMonomialBasis<Displacement<Frame>, Instant,           \
                                       d2, EstrinEvaluator>(                   \
          *dynamic_cast_not_null<                                              \
              PolynomialInMonomialBasis<Displacement<Frame>, Instant,          \
                                        d1, EstrinEvaluator> const*>(          \
                                        polynomial));                          \
    } else {                                                                   \
      LOG(FATAL) << "Inconsistent degrees " << d1 << " and " << d2;            \
    }                                                                          \
  }

template<typename Frame>
template<int aperiodic_degree, int periodic_degree>
PiecewisePoissonSeries<Displacement<Frame>,
                       aperiodic_degree, periodic_degree,
                       EstrinEvaluator>
ContinuousTrajectory<Frame>::ToPiecewisePoissonSeries(
    Instant const& t_min,
    Instant const& t_max) const {
  static_assert(aperiodic_degree >= min_degree &&
                aperiodic_degree <= max_degree);
  // No check on the periodic degree, it plays no role here.
  CHECK(!polynomials_.empty());
  using PiecewisePoisson =
      PiecewisePoissonSeries<Displacement<Frame>,
                             aperiodic_degree, periodic_degree,
                             EstrinEvaluator>;
  using Poisson = PoissonSeries<Displacement<Frame>,
                                aperiodic_degree, periodic_degree,
                                EstrinEvaluator>;

  auto cast_to_degree =
      [](not_null<Polynomial<Displacement<Frame>, Instant> const*> const
             polynomial)
      -> PolynomialInMonomialBasis<Displacement<Frame>, Instant,
                                   aperiodic_degree, EstrinEvaluator> {
    switch (polynomial->degree()) {
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 3, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 4, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 5, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 6, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 7, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 8, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 9, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 10, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 11, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 12, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 13, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 14, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 15, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 16, aperiodic_degree);
      PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS(
          polynomial, 17, aperiodic_degree);
      default:
        LOG(FATAL) << "Unexpected degree " << polynomial->degree();
    }
  };

  std::unique_ptr<PiecewisePoisson> result;

  absl::ReaderMutexLock l(&lock_);
  auto const it_min = FindPolynomialForInstantLocked(t_min);
  auto const it_max = FindPolynomialForInstantLocked(t_max);
  Instant current_t_min = t_min;
  for (auto it = it_min;; ++it) {
    Instant const current_t_max = std::min(t_max, it->t_max);
    Interval<Instant> interval;
    interval.Include(current_t_min);
    interval.Include(current_t_max);
    auto const polynomial_cast_to_degree = cast_to_degree(it->polynomial.get());
    if (result == nullptr) {
      result = std::make_unique<PiecewisePoisson>(
          interval, Poisson(polynomial_cast_to_degree, {{}}));
    } else {
      result->Append(interval, Poisson(polynomial_cast_to_degree, {{}}));
    }
    current_t_min = current_t_max;
    if (it == it_max) {
      break;
    }
  }
  return *result;
}

#undef PRINCIPIA_CAST_TO_POLYNOMIAL_IN_MONOMIAL_BASIS
#endif

template<typename Frame>
void ContinuousTrajectory<Frame>::WriteToMessage(
      not_null<serialization::ContinuousTrajectory*> const message) const {
  absl::ReaderMutexLock l(&lock_);
  CHECK_LT(checkpointer_->oldest_checkpoint(), InfiniteFuture);
  checkpointer_->WriteToMessage(message->mutable_checkpoint());
  step_.WriteToMessage(message->mutable_step());
  tolerance_.WriteToMessage(message->mutable_tolerance());
  polynomial_evaluator_policy_.WriteToMessage(message->mutable_policy());

  // There should be no polynomials before the oldest checkpoint in recent
  // saves, see Ephemeris::AppendMassiveBodiesState.  This has probably been
  // true since Fatou (#2149), but we maintain compatibility with older saves,
  // see #3039.  When such an old save is rewritten, we end up with polynomials
  // before the oldest checkpoint.
  for (auto const& pair : polynomials_) {
    Instant const& t_max = pair.t_max;
    auto const& polynomial = pair.polynomial;
    if (t_max <= checkpointer_->oldest_checkpoint()) {
      auto* const pair = message->add_instant_polynomial_pair();
      t_max.WriteToMessage(pair->mutable_t_max());
      polynomial->WriteToMessage(pair->mutable_polynomial());
    } else {
      break;
    }
  }
  if (first_time_) {
    first_time_->WriteToMessage(message->mutable_first_time());
  }
}

template<typename Frame>
not_null<std::unique_ptr<ContinuousTrajectory<Frame>>>
ContinuousTrajectory<Frame>::ReadFromMessage(
    Instant const& desired_t_min,
    serialization::ContinuousTrajectory const& message)
  requires serializable<Frame> {
  bool const is_pre_cohen = message.series_size() > 0;
  bool const is_pre_fatou = !message.has_checkpoint_time();
  bool const is_pre_grassmann = message.has_adjusted_tolerance() &&
                                message.has_is_unstable() &&
                                message.has_degree() &&
                                message.has_degree_age();
  bool const is_pre_gröbner =
      is_pre_grassmann ||
      (message.instant_polynomial_pair_size() > 0 &&
       message.instant_polynomial_pair(0).polynomial()
           .GetExtension(serialization::PolynomialInMonomialBasis::extension)
           .coefficient(0).has_multivector());
  bool const is_pre_کاشانی = !message.has_policy();
  LOG_IF(WARNING, is_pre_کاشانی)
      << "Reading pre-"
      << (is_pre_cohen       ? "Cohen"
          : is_pre_fatou     ? "Fatou"
          : is_pre_grassmann ? "Grassmann"
          : is_pre_gröbner   ? "Gröbner"
                             : "کاشانی") << " ContinuousTrajectory";

  not_null<std::unique_ptr<ContinuousTrajectory<Frame>>> continuous_trajectory =
      std::make_unique<ContinuousTrajectory<Frame>>(
          Time::ReadFromMessage(message.step()),
          Length::ReadFromMessage(message.tolerance()));
  if (is_pre_cohen) {
    for (auto const& s : message.series()) {
      // Read the polynomial, evaluate it and use the resulting values to build
      // a polynomial in the monomial basis.
      auto const pre_cohen_чебышёв_series =
          PolynomialInЧебышёвBasis<Displacement<Frame>, Instant>::
              ReadFromMessage(s);
      auto const& polynomial = *pre_cohen_чебышёв_series;
      Time const step =
          (polynomial.upper_bound() - polynomial.lower_bound()) / divisions;
      Instant t = polynomial.lower_bound();
      std::vector<Position<Frame>> q;
      std::vector<Velocity<Frame>> v;
      for (int i = 0; i <= divisions; t += step, ++i) {
        q.push_back(polynomial(t) + Frame::origin);
        v.push_back(polynomial.EvaluateDerivative(t));
      }
      Displacement<Frame> error_estimate;  // Should we do something with this?
      continuous_trajectory->polynomials_.emplace_back(
          polynomial.upper_bound(),
          continuous_trajectory->NewhallApproximationInMonomialBasis(
              polynomial.degree(),
              q,
              v,
              polynomial.lower_bound(),
              polynomial.upper_bound(),
              error_estimate));
    }
  } else {
    for (auto const& pair : message.instant_polynomial_pair()) {
      if (is_pre_gröbner) {
        // The easiest way to implement compatibility is to patch the serialized
        // form.
        serialization::Polynomial polynomial = pair.polynomial();
        auto const coefficient0_multivector = polynomial.GetExtension(
            serialization::PolynomialInMonomialBasis::extension).
            coefficient(0).multivector();
        auto* const coefficient0_point = polynomial.MutableExtension(
            serialization::PolynomialInMonomialBasis::extension)->
            mutable_coefficient(0)->mutable_point();
        *coefficient0_point->mutable_multivector() = coefficient0_multivector;

        // Note the use of `EstrinWithoutFMA` below.  FMA was introduced in
        // Grossmann, but unfortunately we didn't notice that it was affecting
        // existing flight plans, and we have no way to tell that a save is
        // exactly pre-Grossmann.  On the other hand, we are sure that pre-
        // Gröbner saves didn't have FMA.  So deserialization is technically
        // incorrect for Gröbner saves, as we will use FMA when reading, even
        // though we didn't have FMA when the save was created.
        continuous_trajectory->polynomials_.emplace_back(
            Instant::ReadFromMessage(pair.t_max()),
            Polynomial<Position<Frame>, Instant>::template ReadFromMessage<
                EstrinWithoutFMA>(polynomial));
      } else {
        serialization::Polynomial const& polynomial = pair.polynomial();
        bool const is_pre_καραθεοδωρή =
            !polynomial.HasExtension(
                serialization::PolynomialInMonomialBasis::extension) ||
            !polynomial
                 .GetExtension(
                     serialization::PolynomialInMonomialBasis::extension)
                 .has_evaluator();
        if (is_pre_καραθεοδωρή) {
          continuous_trajectory->polynomials_.emplace_back(
              Instant::ReadFromMessage(pair.t_max()),
              CanUseHardwareFMA
                  ? Polynomial<Position<Frame>, Instant>::
                        template ReadFromMessage<Estrin>(polynomial)
                  : Polynomial<Position<Frame>, Instant>::
                        template ReadFromMessage<EstrinWithoutFMA>(polynomial));
        } else {
          auto rewritten_polynomial = polynomial;
          if (!CanUseHardwareFMA) {
            // If we are on a machine without FMA, turn evaluators that allow
            // FMA into their counterparts without FMA so that plans made on
            // this machine can be read on a more modern machine.
            // Note that until Lánczos, we created and serialized polynomials
            // with Estrin even when FMA was unavailable, so we cannot trust the
            // evaluator kind here.  Downgrading the evaluators also seems like
            // a better way to handle the case where a save gets moved from a
            // modern machine to an ancient one: in that direction, we cannot
            // preserve the behaviour, but by serializing EstrinWithoutFMA at
            // that transition, whatever work is done on the old machine will be
            // faithfully reproduced when it is read again on a new machine.
            auto const evaluator_kind =
                polynomial
                    .GetExtension(
                        serialization::PolynomialInMonomialBasis::extension)
                    .evaluator()
                    .kind();
            auto& rewritten_evaluator =
                *rewritten_polynomial
                     .MutableExtension(
                         serialization::PolynomialInMonomialBasis::extension)
                     ->mutable_evaluator();
            using Evaluator =
                serialization::PolynomialInMonomialBasis::Evaluator;
            switch (evaluator_kind) {
              case Evaluator::ESTRIN:
                rewritten_evaluator.set_kind(Evaluator::ESTRIN_WITHOUT_FMA);
                break;
              case Evaluator::HORNER:
                rewritten_evaluator.set_kind(Evaluator::HORNER_WITHOUT_FMA);
                break;
              case Evaluator::HORNER_WITHOUT_FMA:
              case Evaluator::ESTRIN_WITHOUT_FMA:
                break;
            }
          }
          continuous_trajectory->polynomials_.emplace_back(
              Instant::ReadFromMessage(pair.t_max()),
              Polynomial<Position<Frame>, Instant>::ReadFromMessage(
                  rewritten_polynomial));
        }
      }
    }
  }
  if (is_pre_کاشانی) {
    // See the comment above for the defaults here.
    if (is_pre_gröbner || !CanUseHardwareFMA) {
      continuous_trajectory->polynomial_evaluator_policy_ =
          Policy::AlwaysEstrinWithoutFMA();
    } else {
      continuous_trajectory->polynomial_evaluator_policy_ =
          Policy::AlwaysEstrin();
    }
  } else {
    continuous_trajectory->polynomial_evaluator_policy_ =
        CanUseHardwareFMA ? Policy::ReadFromMessage(message.policy())
                          : Policy::AlwaysEstrinWithoutFMA();
  }
  if (message.has_first_time()) {
    continuous_trajectory->first_time_ =
        Instant::ReadFromMessage(message.first_time());
  }

  if (is_pre_grassmann) {
    serialization::ContinuousTrajectory serialized_continuous_trajectory;
    auto* const checkpoint = serialized_continuous_trajectory.add_checkpoint();
    if (is_pre_fatou) {
      *checkpoint->mutable_time() =
          message.last_point()[message.last_point_size() - 1].instant();
    } else {
      *checkpoint->mutable_time() = message.checkpoint_time();
    }
    *checkpoint->mutable_adjusted_tolerance() = message.adjusted_tolerance();
    checkpoint->set_is_unstable(message.is_unstable());
    checkpoint->set_degree(message.degree());
    checkpoint->set_degree_age(message.degree_age());
    for (const auto& last_point : message.last_point()) {
      *checkpoint->add_last_point() = last_point;
    }
    continuous_trajectory->checkpointer_ =
        Checkpointer<serialization::ContinuousTrajectory>::ReadFromMessage(
            continuous_trajectory->MakeCheckpointerWriter(),
            continuous_trajectory->MakeCheckpointerReader(),
            serialized_continuous_trajectory.checkpoint());
  } else {
    continuous_trajectory->checkpointer_ =
        Checkpointer<serialization::ContinuousTrajectory>::ReadFromMessage(
            continuous_trajectory->MakeCheckpointerWriter(),
            continuous_trajectory->MakeCheckpointerReader(),
            message.checkpoint());
  }

  // There should always be a checkpoint, either at the end of the trajectory,
  // in the pre-Grassmann compatibility case; or at the first point of the
  // trajectory, for modern saves (see the comment in WriteToMessage).  In the
  // pre-Grassman case the (only) checkpoint may be after `desired_t_min`, but
  // then we should have polynomials covering that time, see #3039.
  auto const status =
      continuous_trajectory->checkpointer_->ReadFromCheckpointAtOrBefore(
          desired_t_min);
  if (!status.ok()) {
    CHECK(absl::IsNotFound(status));
    CHECK_LE(continuous_trajectory->t_min(), desired_t_min);
  }

  return continuous_trajectory;
}

template<typename Frame>
void ContinuousTrajectory<Frame>::WriteToCheckpoint(Instant const& t) const {
  checkpointer_->WriteToCheckpoint(t);
}

template<typename Frame>
absl::Status ContinuousTrajectory<Frame>::ReadFromCheckpointAt(
    Instant const& t,
    Checkpointer<serialization::ContinuousTrajectory>::Reader const& reader)
    const {
  return checkpointer_->ReadFromCheckpointAt(t, reader);
}

template<typename Frame>
Checkpointer<serialization::ContinuousTrajectory>::Writer
ContinuousTrajectory<Frame>::MakeCheckpointerWriter() {
  if constexpr (serializable<Frame>) {
    return [this](
        not_null<
            serialization::ContinuousTrajectory::Checkpoint*> const message) {
      absl::ReaderMutexLock l(&lock_);
      adjusted_tolerance_.WriteToMessage(message->mutable_adjusted_tolerance());
      message->set_is_unstable(is_unstable_);
      message->set_degree(degree_);
      message->set_degree_age(degree_age_);
      for (auto const& [instant, degrees_of_freedom] : last_points_) {
        not_null<serialization::ContinuousTrajectory::
                     InstantaneousDegreesOfFreedom*> const
            instantaneous_degrees_of_freedom = message->add_last_point();
        instant.WriteToMessage(
            instantaneous_degrees_of_freedom->mutable_instant());
        degrees_of_freedom.WriteToMessage(
            instantaneous_degrees_of_freedom->mutable_degrees_of_freedom());
      }
    };
  } else {
    return nullptr;
  }
}

template<typename Frame>
Checkpointer<serialization::ContinuousTrajectory>::Reader
ContinuousTrajectory<Frame>::MakeCheckpointerReader() {
  if constexpr (serializable<Frame>) {
    return [this](
               serialization::ContinuousTrajectory::Checkpoint const& message) {
      absl::MutexLock l(&lock_);

      // Load the members that are recorded in the checkpoint.
      adjusted_tolerance_ =
          Length::ReadFromMessage(message.adjusted_tolerance());
      is_unstable_ = message.is_unstable();
      degree_ = message.degree();
      degree_age_ = message.degree_age();
      last_points_.clear();
      for (auto const& l : message.last_point()) {
        last_points_.push_back(
            {Instant::ReadFromMessage(l.instant()),
             DegreesOfFreedom<Frame>::ReadFromMessage(l.degrees_of_freedom())});
      }

      // Restore the other members to their state at the time of the checkpoint.
      if (last_points_.empty()) {
        polynomials_.clear();
        first_time_ = std::nullopt;
      } else {
        // Locate the polynomial that ends at the first last_point_.  Note that
        // we cannot use FindPolynomialForInstantLocked because it calls
        // lower_bound and we don't want to change its behaviour.
        Instant const& oldest_time = last_points_.front().first;
        // If oldest_time is the t_max of some polynomial, then the returned
        // iterator points to the next polynomial.
        auto const it =
            std::upper_bound(polynomials_.begin(),
                             polynomials_.end(),
                             oldest_time,
                             [](Instant const& left,
                                InstantPolynomialPair const& right) {
                               return left < right.t_max;
                             });
        polynomials_.erase(it, polynomials_.end());
        if (polynomials_.empty()) {
          first_time_ = oldest_time;
        }
      }
      last_accessed_polynomial_ = 0;  // Always a valid value.

      return absl::OkStatus();
    };
  } else {
    return nullptr;
  }
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_min_locked() const {
  if (polynomials_.empty()) {
    return InfiniteFuture;
  }
  return *first_time_;
}

template<typename Frame>
Instant ContinuousTrajectory<Frame>::t_max_locked() const {
  if (polynomials_.empty()) {
    return InfinitePast;
  }
  return polynomials_.crbegin()->t_max;
}

template<typename Frame>
Position<Frame> ContinuousTrajectory<Frame>::EvaluatePositionLocked(
    Instant const& time) const {
  CHECK_LE(t_min_locked(), time);
  CHECK_GE(t_max_locked(), time);
  auto const it = FindPolynomialForInstantLocked(time);
  CHECK(it != polynomials_.end());
  auto const& polynomial = *it->polynomial;
  return polynomial(time);
}

template<typename Frame>
Velocity<Frame> ContinuousTrajectory<Frame>::EvaluateVelocityLocked(
    Instant const& time) const {
  CHECK_LE(t_min_locked(), time);
  CHECK_GE(t_max_locked(), time);
  auto const it = FindPolynomialForInstantLocked(time);
  CHECK(it != polynomials_.end());
  auto const& polynomial = *it->polynomial;
  return polynomial.EvaluateDerivative(time);
}

template<typename Frame>
DegreesOfFreedom<Frame>
ContinuousTrajectory<Frame>::EvaluateDegreesOfFreedomLocked(
    Instant const& time) const {
  CHECK_LE(t_min_locked(), time);
  CHECK_GE(t_max_locked(), time);
  auto const it = FindPolynomialForInstantLocked(time);
  CHECK(it != polynomials_.end());
  auto const& polynomial = *it->polynomial;
  return DegreesOfFreedom<Frame>(polynomial(time),
                                 polynomial.EvaluateDerivative(time));
}

template<typename Frame>
ContinuousTrajectory<Frame>::ContinuousTrajectory()
    : checkpointer_(
          make_not_null_unique<
              Checkpointer<serialization::ContinuousTrajectory>>(
          /*reader=*/nullptr,
          /*writer=*/nullptr)),
      polynomial_evaluator_policy_(CanUseHardwareFMA
                                       ? Policy::AlwaysEstrin()
                                       : Policy::AlwaysEstrinWithoutFMA()) {}

template<typename Frame>
ContinuousTrajectory<Frame>::InstantPolynomialPair::InstantPolynomialPair(
    Instant const t_max,
    not_null<std::unique_ptr<Polynomial<Position<Frame>, Instant>>>
        polynomial)
    : t_max(t_max),
      polynomial(std::move(polynomial)) {}

template<typename Frame>
not_null<std::unique_ptr<Polynomial<Position<Frame>, Instant>>>
ContinuousTrajectory<Frame>::NewhallApproximationInMonomialBasis(
    int degree,
    std::vector<Position<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v,
    Instant const& t_min,
    Instant const& t_max,
    Displacement<Frame>& error_estimate) const {
  return numerics::_newhall::NewhallApproximationInMonomialBasis<
             Position<Frame>>(degree,
                              q, v,
                              t_min, t_max,
                              polynomial_evaluator_policy_,
                              error_estimate);
}

template<typename Frame>
absl::Status ContinuousTrajectory<Frame>::ComputeBestNewhallApproximation(
    Instant const& time,
    std::vector<Position<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v) {
  lock_.AssertHeld();
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

  // Estimate the error.  For initializing `previous_error_estimate`, any value
  // greater than `error_estimate` will do.
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
    return absl::OkStatus();
  } else {
    std::stringstream message;
    message << "Error trying to fit a smooth polynomial to the trajectory. "
            << "The approximation error jumped from "
            << previous_adjusted_tolerance << " to " << adjusted_tolerance_
            << " at time " << time << ". The last position is " << q.back()
            << " and the last velocity is " << v.back()
            << ". An apocalypse occurred and two celestials probably "
            << "collided because your solar system is unstable.";
    return absl::InvalidArgumentError(message.str());
  }
}

template<typename Frame>
typename ContinuousTrajectory<Frame>::InstantPolynomialPairs::const_iterator
ContinuousTrajectory<Frame>::FindPolynomialForInstantLocked(
    Instant const& time) const {
  // This returns the first polynomial `p` such that `time <= p.t_max`.
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

}  // namespace internal
}  // namespace _continuous_trajectory
}  // namespace physics
}  // namespace principia
