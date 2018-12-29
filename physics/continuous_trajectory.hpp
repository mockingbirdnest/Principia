﻿
#pragma once

#include <atomic>
#include <optional>
#include <utility>
#include <vector>

#include "base/arena_allocator.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/polynomial.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_continuous_trajectory {

using base::ArenaAllocator;
using base::not_null;
using base::Status;
using geometry::Displacement;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using quantities::Length;
using quantities::Time;
using numerics::Polynomial;

template<typename Frame>
class TestableContinuousTrajectory;

template<typename Frame>
class ContinuousTrajectory : public Trajectory<Frame> {
 public:
  // A |Checkpoint| contains the impermanent state of a trajectory, i.e., the
  // state that gets incrementally updated as the polynomials are constructed.
  // The client may get a |Checkpoint| at any time and use it to serialize the
  // trajectory up to and including the time designated by the |Checkpoint|.
  class Checkpoint;

  // Constructs a trajectory with the given time |step|.  Because the Чебышёв
  // polynomials have values in the range [-1, 1], the error resulting of
  // truncating the infinite Чебышёв series to a finite degree are a small
  // multiple of the coefficient of highest degree (assuming that the series
  // converges reasonably well).  Thus, we pick the degree of the series so that
  // the coefficient of highest degree is less than |tolerance|.
  ContinuousTrajectory(Time const& step,
                       Length const& tolerance);
  virtual ~ContinuousTrajectory() = default;

  ContinuousTrajectory(ContinuousTrajectory const&) = delete;
  ContinuousTrajectory(ContinuousTrajectory&&) = delete;
  ContinuousTrajectory& operator=(ContinuousTrajectory const&) = delete;
  ContinuousTrajectory& operator=(ContinuousTrajectory&&) = delete;

  // Returns true iff this trajectory cannot be evaluated for any time.
  bool empty() const;

  // The average degree of the polynomials for the trajectory.  Only useful for
  // benchmarking or analyzing performance.  Do not use in real code.
  double average_degree() const;

  // Appends one point to the trajectory.  |time| must be after the last time
  // passed to |Append| if the trajectory is not empty.  The |time|s passed to
  // successive calls to |Append| must be equally spaced with the |step| given
  // at construction.
  Status Append(Instant const& time,
                DegreesOfFreedom<Frame> const& degrees_of_freedom);

  // Removes all data for times strictly less than |time|.
  void ForgetBefore(Instant const& time);

  // Implementation of the interface |Trajectory|.

  // |t_max| may be less than the last time passed to Append.  For an empty
  // trajectory, an infinity with the proper sign is returned.
  Instant t_min() const override;
  Instant t_max() const override;

  Position<Frame> EvaluatePosition(Instant const& time) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& time) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& time) const override;

  // End of the implementation of the interface.

  // Returns a checkpoint for the current state of this object.
  Checkpoint GetCheckpoint() const;

  // Serializes the current state of this object.
  void WriteToMessage(
      not_null<serialization::ContinuousTrajectory*> message) const;
  // Serializes the state of this object as it existed when the checkpoint was
  // taken.
  void WriteToMessage(not_null<serialization::ContinuousTrajectory*> message,
                      Checkpoint const& checkpoint) const;
  static not_null<std::unique_ptr<ContinuousTrajectory>> ReadFromMessage(
      serialization::ContinuousTrajectory const& message);

  // A |Checkpoint| contains the impermanent state of a trajectory, i.e., the
  // state that gets incrementally updated as the polynomials are constructed.
  // The client may get a |Checkpoint| at any time and use it to serialize the
  // trajectory up to and including the time designated by the |Checkpoint|.
  // The only thing that clients may do with |Checkpoint| objects is to
  // initialize them with GetCheckpoint.
  class Checkpoint final {
   public:
    // Returns true if this checkpoint is after |time| and would remain valid
    // after a call to |ForgetBefore(time)|.
    bool IsAfter(Instant const& time) const;

   private:
    // The members have the same meaning as those of class
    // |ContinuousTrajectory|.
    Checkpoint(Instant const& t_max,
               Length const& adjusted_tolerance,
               bool is_unstable,
               int degree,
               int degree_age,
               std::vector<std::pair<Instant, DegreesOfFreedom<Frame>>> const&
                   last_points);
    Instant t_max_;
    Length adjusted_tolerance_;
    bool is_unstable_;
    int degree_;
    int degree_age_;
    std::vector<std::pair<Instant, DegreesOfFreedom<Frame>>> last_points_;
    friend class ContinuousTrajectory<Frame>;
  };

 protected:
  // For mocking.
  ContinuousTrajectory();

 private:
  // Each polynomial is valid over an interval [t_min, t_max].  Polynomials are
  // stored in this vector sorted by their |t_max|, as it turns out that we
  // never need to extract their |t_min|.  Logically, the |t_min| for a
  // polynomial is the |t_max| of the previous one.  The first polynomial has a
  // |t_min| which is |*first_time_|.  The pointer is owned by |arena_|;
  struct InstantPolynomialPair {
    InstantPolynomialPair(
        Instant t_max,
        not_null<Polynomial<Displacement<Frame>, Instant>*> polynomial);
    Instant t_max;
    not_null<Polynomial<Displacement<Frame>, Instant>*> polynomial;
  };
  using InstantPolynomialPairs =
      std::vector<InstantPolynomialPair, ArenaAllocator<InstantPolynomialPair>>;

  // May be overridden for testing.
  virtual not_null<Polynomial<Displacement<Frame>, Instant>*>
  NewhallApproximationInMonomialBasis(
      int degree,
      std::vector<Displacement<Frame>> const& q,
      std::vector<Velocity<Frame>> const& v,
      Instant const& t_min,
      Instant const& t_max,
      Displacement<Frame>& error_estimate);

  // Computes the best Newhall approximation based on the desired tolerance.
  // Adjust the |degree_| and other member variables to stay within the
  // tolerance while minimizing the computational cost and avoiding numerical
  // instabilities.
  Status ComputeBestNewhallApproximation(
      Instant const& time,
      std::vector<Displacement<Frame>> const& q,
      std::vector<Velocity<Frame>> const& v);

  // Returns an iterator to the polynomial applicable for the given |time|, or
  // |begin()| if |time| is before the first polynomial or |end()| if |time| is
  // after the last polynomial.  Time complexity is O(N Log N).
  typename InstantPolynomialPairs::const_iterator
  FindPolynomialForInstant(Instant const& time) const;

  // Construction parameters;
  Time const step_;
  Length const tolerance_;

  // Initially set to the construction parameters, and then adjusted when we
  // choose the degree.
  Length adjusted_tolerance_;
  bool is_unstable_;

  // The degree of the approximation and its age in number of Newhall
  // approximations.
  int degree_;
  int degree_age_;

  // Storage management.
  google::protobuf::Arena arena_;
  ArenaAllocator<InstantPolynomialPairs> allocator_;

  // The polynomials are in increasing time order.
  InstantPolynomialPairs polynomials_;

  // Lookups into |polynomials_| are expensive because they entail a binary
  // search into a vector that grows over time.  In benchmarks, this can be as
  // costly as the polynomial evaluation itself.  The accesses are not random,
  // though, they are clustered in time and (slowly) increasing.  To take
  // advantage of this, we keep track of the index of the last accessed
  // polynomial and first try to see if the new lookup is for the same
  // polynomial.  This makes us O(1) instead of O(Log N) most of the time and it
  // speeds up the lookup by a factor of 7.  This member is atomic because of
  // multithreading in the ephemeris, and is mutable to maintain the fiction
  // that evaluation has no side effects.  In the presence of multithreading it
  // may be that different threads would want to access polynomials at different
  // indices, but by and large the threads progress in parallel, and benchmarks
  // show that there is no adverse performance effects.  Any value in the range
  // of |polynomials_| or 0 is correct.
  mutable std::atomic_int last_accessed_polynomial_ = 0;

  // The time at which this trajectory starts.  Set for a nonempty trajectory.
  std::optional<Instant> first_time_;

  // The points that have not yet been incorporated in a polynomial.  Nonempty
  // for a nonempty trajectory.
  // |last_points_.begin()->first == polynomials_.back().t_max|
  std::vector<std::pair<Instant, DegreesOfFreedom<Frame>>> last_points_;

  friend class TestableContinuousTrajectory<Frame>;
};

}  // namespace internal_continuous_trajectory

using internal_continuous_trajectory::ContinuousTrajectory;

}  // namespace physics
}  // namespace principia

#if !PHYSICS_DLL_IMPORT
#include "physics/continuous_trajectory_body.hpp"
#endif
