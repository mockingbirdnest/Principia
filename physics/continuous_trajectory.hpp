#pragma once

#include <vector>
#include <utility>

#include "geometry/named_quantities.hpp"
#include "numerics/чебышёв_series.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {

using geometry::Instant;
using quantities::Length;
using quantities::Time;
using numerics::ЧебышёвSeries;

namespace physics {

template<typename Frame>
class ContinuousTrajectory {
 public:
  // A |Hint| is used to speed up the evaluation of trajectories.  When
  // repeatedly calling one of the evaluation functions with increasing values
  // of the |time| parameter, evaluation may be faster if the same |Hint| object
  // is passed to all the calls.
  class Hint;

  // Constructs a trajectory with the given time |step|.  Because the Чебышёв
  // polynomials have values in the range [-1, 1], the error resulting of
  // truncating the infinite Чебышёв series to a finite degree are a small
  // multiple of the coefficient of highest degree (assuming that the series
  // converges reasonably well).  Thus, we pick the degree of the series so that
  // the coefficient of highest degree is less than |tolerance|.
  ContinuousTrajectory(Time const& step,
                       Length const& tolerance);
  ~ContinuousTrajectory() = default;

  ContinuousTrajectory(ContinuousTrajectory const&) = delete;
  ContinuousTrajectory(ContinuousTrajectory&&) = delete;
  ContinuousTrajectory& operator=(ContinuousTrajectory const&) = delete;
  ContinuousTrajectory& operator=(ContinuousTrajectory&&) = delete;

  // Returns true iff this trajectory cannot be evaluated for any time.
  bool empty() const;

  // The time range for which the trajectory can be evaluated.  Note that
  // |t_max| may be less than the last time passed to Append.  For an empty
  // trajectory, an infinity with the proper sign is returned.
  Instant t_min() const;
  Instant t_max() const;

  // Appends one point to the trajectory.  |time| must be after the last time
  // passed to |Append| if the trajectory is not empty.  The |time|s passed to
  // successive calls to |Append| must be equally spaced with the |step| given
  // at construction.
  void Append(Instant const& time,
              DegreesOfFreedom<Frame> const& degrees_of_freedom);

  // Removes all data for times strictly less than |time|.
  void ForgetBefore(Instant const& time);

  // Evaluates the trajectory at the given |time|, which must be in
  // [t_min(), t_max()].  The |hint| may be used to speed up evaluation
  // in increasing time order.  It may be a nullptr (in which case no speed-up
  // takes place).
  Position<Frame> EvaluatePosition(Instant const& time,
                                   Hint* const hint) const;
  Velocity<Frame> EvaluateVelocity(Instant const& time,
                                   Hint* const hint) const;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(Instant const& time,
                                                   Hint* const hint) const;

  void WriteToMessage(
      not_null<serialization::ContinuousTrajectory*> const message) const;
  static not_null<std::unique_ptr<ContinuousTrajectory>> ReadFromMessage(
      serialization::ContinuousTrajectory const& message);

  // The only thing that clients may do with |Hint| objects is to
  // default-initialize them.
  class Hint {
   public:
    Hint();
   private:
    int index_;
    friend class ContinuousTrajectory<Frame>;
  };

 private:
  // Computes the best Newhall approximation based on the desired tolerance.
  // Adjust the |degree_| and other member variables to stay within the
  // tolerance while minimizing the computational cost and avoiding numerical
  // instabilities.
  void ComputeBestNewhallApproximation(
      Instant const& time,
      std::vector<Displacement<Frame>> const& q,
      std::vector<Velocity<Frame>> const& v,
      ЧебышёвSeries<Displacement<Frame>> (*newhall_approximation)(
          int const degree,
          std::vector<Displacement<Frame>> const& q,
          std::vector<Velocity<Frame>> const& v,
          Instant const& t_min,
          Instant const& t_max));

  // Returns an iterator to the series applicable for the given |time|, or
  // |begin()| if |time| is before the first series or |end()| if |time| is
  // after the last series.  Time complexity is O(N Log N).
  typename std::vector<ЧебышёвSeries<Displacement<Frame>>>::const_iterator
  FindSeriesForInstant(Instant const& time) const;

  // Returns true if the given |hint| is usable for the given |time|.  If it is,
  // |hint->index| is the index of the series to use.
  bool MayUseHint(Instant const& time, Hint* const hint) const;

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

  // The series are in increasing time order.  Their intervals are consecutive.
  std::vector<ЧебышёвSeries<Displacement<Frame>>> series_;

  // The time at which this trajectory starts.  Set for a nonempty trajectory.
  // |first_time_ >= series_.front().t_min()|
  std::unique_ptr<Instant> first_time_;  // std::optional.

  // The points that have not yet been incorporated in a series.  Nonempty for a
  // nonempty trajectory.
  // |last_points_.begin()->first == series_.back().t_max()|
  std::vector<std::pair<Instant, DegreesOfFreedom<Frame>>> last_points_;

  friend class ContinuousTrajectoryTest;
  friend class TransformsTest;
};

}  // namespace physics
}  // namespace principia

#include "physics/continuous_trajectory_body.hpp"
