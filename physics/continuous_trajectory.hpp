
#pragma once

#include <atomic>
#include <optional>
#include <utility>
#include <vector>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/poisson_series.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "physics/checkpointer.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_continuous_trajectory {

using base::not_null;
using base::Status;
using geometry::Displacement;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using quantities::Length;
using quantities::Time;
using numerics::EstrinEvaluator;
using numerics::PiecewisePoissonSeries;
using numerics::Polynomial;

template<typename Frame>
class TestableContinuousTrajectory;

// This class is thread-safe, but the client must be aware that if, for
// instance, the trajectory is appended to asynchronously, successive calls to
// |t_max()| may return different values.
template<typename Frame>
class ContinuousTrajectory : public Trajectory<Frame> {
 public:
  // Constructs a trajectory with the given time |step|.  Because the Чебышёв
  // polynomials have values in the range [-1, 1], the error resulting of
  // truncating the infinite Чебышёв series to a finite degree are a small
  // multiple of the coefficient of highest degree (assuming that the series
  // converges reasonably well).  Thus, we pick the degree of the series so that
  // the coefficient of highest degree is less than |tolerance|.
  ContinuousTrajectory(Time const& step,
                       Length const& tolerance);

  ContinuousTrajectory(ContinuousTrajectory const&) = delete;
  ContinuousTrajectory(ContinuousTrajectory&&) = delete;
  ContinuousTrajectory& operator=(ContinuousTrajectory const&) = delete;
  ContinuousTrajectory& operator=(ContinuousTrajectory&&) = delete;

  // Returns true iff this trajectory cannot be evaluated for any time.
  bool empty() const EXCLUDES(lock_);

  // The average degree of the polynomials for the trajectory.  Only useful for
  // benchmarking or analyzing performance.  Do not use in real code.
  double average_degree() const EXCLUDES(lock_);

  // Appends one point to the trajectory.  |time| must be after the last time
  // passed to |Append| if the trajectory is not empty.  The |time|s passed to
  // successive calls to |Append| must be equally spaced with the |step| given
  // at construction.
  Status Append(Instant const& time,
                DegreesOfFreedom<Frame> const& degrees_of_freedom)
      EXCLUDES(lock_);

  // Removes all data for times strictly less than |time|.
  void ForgetBefore(Instant const& time) EXCLUDES(lock_);

  // Implementation of the interface |Trajectory|.

  // |t_max| may be less than the last time passed to Append.  For an empty
  // trajectory, an infinity with the proper sign is returned.
  Instant t_min() const override EXCLUDES(lock_);
  Instant t_max() const override EXCLUDES(lock_);

  Position<Frame> EvaluatePosition(Instant const& time) const override
      EXCLUDES(lock_);
  Velocity<Frame> EvaluateVelocity(Instant const& time) const override
      EXCLUDES(lock_);
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& time) const override EXCLUDES(lock_);

  // End of the implementation of the interface.

  // Returns the degree for a piecewise Poisson series covering the given time
  // interval.
  int PiecewisePoissonSeriesDegree(Instant const& t_min,
                                   Instant const& t_max) const;

  // Computes a piecewise Poisson series covering the given time interval.  The
  // degree must be at least the one returned by the preceding function.
  template<int degree>
  PiecewisePoissonSeries<Displacement<Frame>, degree, EstrinEvaluator>
  ToPiecewisePoissonSeries(Instant const& t_min,
                           Instant const& t_max,
                           Instant const& origin) const;

  void WriteToMessage(not_null<serialization::ContinuousTrajectory*> message)
      const EXCLUDES(lock_);
  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  static not_null<std::unique_ptr<ContinuousTrajectory>> ReadFromMessage(
      serialization::ContinuousTrajectory const& message);

  // Checkpointing support.  The checkpointer is exposed to make it possible for
  // Ephemeris to create synchronized checkpoints of its state and that of its
  // trajectories.
  Checkpointer<serialization::ContinuousTrajectory>& checkpointer();
  void WriteToCheckpoint(
      not_null<serialization::ContinuousTrajectory*> message);
  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  bool ReadFromCheckpoint(serialization::ContinuousTrajectory const& message);

 protected:
  // For mocking.
  ContinuousTrajectory();

 private:
  // Each polynomial is valid over an interval [t_min, t_max].  Polynomials are
  // stored in this vector sorted by their |t_max|, as it turns out that we
  // never need to extract their |t_min|.  Logically, the |t_min| for a
  // polynomial is the |t_max| of the previous one.  The first polynomial has a
  // |t_min| which is |*first_time_|.
  // TODO(phl): These should be polynomials returning Position<Frame>.
  struct InstantPolynomialPair {
    InstantPolynomialPair(
        Instant t_max,
        not_null<std::unique_ptr<Polynomial<Displacement<Frame>, Instant>>>
            polynomial);
    Instant t_max;
    not_null<std::unique_ptr<Polynomial<Displacement<Frame>, Instant>>>
        polynomial;
  };
  using InstantPolynomialPairs = std::vector<InstantPolynomialPair>;

  Instant t_min_locked() const REQUIRES_SHARED(lock_);
  Instant t_max_locked() const REQUIRES_SHARED(lock_);

  // Really a static method, but may be overridden for testing.
  virtual not_null<std::unique_ptr<Polynomial<Displacement<Frame>, Instant>>>
  NewhallApproximationInMonomialBasis(
      int degree,
      std::vector<Displacement<Frame>> const& q,
      std::vector<Velocity<Frame>> const& v,
      Instant const& t_min,
      Instant const& t_max,
      Displacement<Frame>& error_estimate) const;

  // Computes the best Newhall approximation based on the desired tolerance.
  // Adjust the |degree_| and other member variables to stay within the
  // tolerance while minimizing the computational cost and avoiding numerical
  // instabilities.
  Status ComputeBestNewhallApproximation(
      Instant const& time,
      std::vector<Displacement<Frame>> const& q,
      std::vector<Velocity<Frame>> const& v) REQUIRES(lock_);

  // Returns an iterator to the polynomial applicable for the given |time|, or
  // |begin()| if |time| is before the first polynomial or |end()| if |time| is
  // after the last polynomial.  Time complexity is O(N Log N).
  typename InstantPolynomialPairs::const_iterator
  FindPolynomialForInstant(Instant const& time) const REQUIRES_SHARED(lock_);

  // Construction parameters;
  Time const step_;
  Length const tolerance_;
  Checkpointer<serialization::ContinuousTrajectory> checkpointer_;

  mutable absl::Mutex lock_;

  // Initially set to the construction parameters, and then adjusted when we
  // choose the degree.
  Length adjusted_tolerance_ GUARDED_BY(lock_);
  bool is_unstable_ GUARDED_BY(lock_);

  // The degree of the approximation and its age in number of Newhall
  // approximations.
  int degree_ GUARDED_BY(lock_);
  int degree_age_ GUARDED_BY(lock_);

  // The polynomials are in increasing time order.
  InstantPolynomialPairs polynomials_ GUARDED_BY(lock_);

  // Lookups into |polynomials_| are expensive because they entail a binary
  // search into a vector that grows over time.  In benchmarks, this can be as
  // costly as the polynomial evaluation itself.  The accesses are not random,
  // though, they are clustered in time and (slowly) increasing.  To take
  // advantage of this, we keep track of the index of the last accessed
  // polynomial and first try to see if the new lookup is for the same
  // polynomial.  This makes us O(1) instead of O(Log N) most of the time and it
  // speeds up the lookup by a factor of 7.  This member is mutable to maintain
  // the fiction that evaluation has no side effects.  In the presence of
  // multithreading it may be that different threads would want to access
  // polynomials at different indices, but by and large the threads progress in
  // parallel, and benchmarks show that there is no adverse performance effects.
  // Any value in the range of |polynomials_| or 0 is correct.
  mutable std::int64_t last_accessed_polynomial_ GUARDED_BY(lock_) = 0;

  // The time at which this trajectory starts.  Set for a nonempty trajectory.
  std::optional<Instant> first_time_ GUARDED_BY(lock_);

  // The points that have not yet been incorporated in a polynomial.  Nonempty
  // for a nonempty trajectory.
  // |last_points_.begin()->first == polynomials_.back().t_max|
  std::vector<std::pair<Instant, DegreesOfFreedom<Frame>>> last_points_
      GUARDED_BY(lock_);

  friend class TestableContinuousTrajectory<Frame>;
};

}  // namespace internal_continuous_trajectory

using internal_continuous_trajectory::ContinuousTrajectory;

}  // namespace physics
}  // namespace principia

#include "physics/continuous_trajectory_body.hpp"
