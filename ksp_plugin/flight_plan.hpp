
#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace internal_flight_plan {

using base::Error;
using base::not_null;
using base::Status;
using geometry::Instant;
using integrators::AdaptiveStepSizeIntegrator;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using quantities::Length;
using quantities::Mass;
using quantities::Speed;

// A chain of trajectories obtained by executing the corresponding
// |NavigationManœuvre|s.
class FlightPlan {
 public:
  // Creates a |FlightPlan| with no burns starting at |initial_time| with
  // |initial_degrees_of_freedom| and with the given |initial_mass|.  The
  // trajectories are computed using the given parameters by the given
  // |ephemeris|.  The flight plan contains a single coast which, if possible
  // ends at |desired_final_time|.
  FlightPlan(Mass const& initial_mass,
             Instant const& initial_time,
             DegreesOfFreedom<Barycentric> const& initial_degrees_of_freedom,
             Instant const& desired_final_time,
             not_null<Ephemeris<Barycentric>*> ephemeris,
             Ephemeris<Barycentric>::AdaptiveStepParameters const&
                 adaptive_step_parameters,
             Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
                 generalized_adaptive_step_parameters);
  virtual ~FlightPlan() = default;

  // Construction parameters.
  virtual Instant initial_time() const;
  virtual Instant desired_final_time() const;

  // End time of the last coast.  If this is less than |desired_final_time()|,
  // there is at least an anomalous manœuvre.
  virtual Instant actual_final_time() const;

  // The number of manœuvres in the flight plan.
  virtual int number_of_manœuvres() const;

  // The number of manœuvres at the end of flight plan that are anomalous, i.e.,
  // lead to their burn or a subsequent trajectory failing to integrate.  Note
  // that this returns 1 for a flight plan without manœuvres if the first coast
  // fails to integrate.  The functions that change manœuvres may change the
  // number of anomalous manœuvres.
  virtual int number_of_anomalous_manœuvres() const;

  // Returns the specified manœuvre.  |index| must be in
  // [0, number_of_manœuvres()[.
  virtual NavigationManœuvre const& GetManœuvre(int index) const;

  // Appends a manœuvre using the specified |burn|.  |index| must be in
  // [0, number_of_manœuvres()[.  Returns an error and has no effect if the
  // given |burn| cannot fit between the preceding burn and the end of the
  // flight plan.  Otherwise, updates the flight plan and returns the
  // integration status.
  virtual Status Append(NavigationManœuvre::Burn const& burn);

  // Forgets the flight plan at least before |time|.  The actual cutoff time
  // will be in a coast trajectory and may be after |time|.  |on_empty| is run
  // if the flight plan would become empty (it is not modified before running
  // |on_empty|).
  virtual void ForgetBefore(Instant const& time,
                            std::function<void()> const& on_empty);


  // Removes the last manœuvre.
  virtual Status RemoveLast();

  // Replaces a manœuvre with one using the specified |burn|.  |index| must be
  // in [0, number_of_manœuvres()[.  Returns an error and has no effect if the
  // given |burn| cannot fit between the preceding and following burns.
  // Otherwise, updates the flight plan and returns the integration status.
  virtual Status Replace(NavigationManœuvre::Burn const& burn, int index);

  // Same as above, but for the last manœuvre.  |number_of_manœuvres()| must be
  // at least one..
  virtual Status ReplaceLast(NavigationManœuvre::Burn const& burn);

  // Updates the desired final time of the flight plan.  Returns an error and
  // has no effect |desired_final_time| is before the beginning of the last
  // coast.
  virtual Status SetDesiredFinalTime(Instant const& desired_final_time);

  // Sets the parameters used to compute the trajectories and recomputes them
  // all.  Returns the integration status.
  virtual Status SetAdaptiveStepParameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          adaptive_step_parameters,
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
          generalized_adaptive_step_parameters);

  virtual Ephemeris<Barycentric>::AdaptiveStepParameters const&
  adaptive_step_parameters() const;
  virtual Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
  generalized_adaptive_step_parameters() const;

  // Returns the number of trajectories in this object.
  virtual int number_of_segments() const;

  // |index| must be in [0, number_of_segments()[.  Sets the iterators to denote
  // the given trajectory.
  virtual void GetSegment(int index,
                          DiscreteTrajectory<Barycentric>::Iterator& begin,
                          DiscreteTrajectory<Barycentric>::Iterator& end) const;
  virtual void GetAllSegments(
      DiscreteTrajectory<Barycentric>::Iterator& begin,
      DiscreteTrajectory<Barycentric>::Iterator& end) const;

  void WriteToMessage(not_null<serialization::FlightPlan*> message) const;

  // This may return a null pointer if the flight plan contained in the
  // |message| is anomalous.
  static std::unique_ptr<FlightPlan> ReadFromMessage(
      serialization::FlightPlan const& message,
      not_null<Ephemeris<Barycentric>*> ephemeris);

  static constexpr std::int64_t max_ephemeris_steps_per_frame = 1000;

  static constexpr Error bad_desired_final_time = Error::INVALID_ARGUMENT;
  static constexpr Error does_not_fit = Error::INVALID_ARGUMENT;
  static constexpr Error singular = Error::INVALID_ARGUMENT;

 protected:
  // For mocking.
  FlightPlan();

 private:
  // Clears and recomputes all trajectories in |segments_|.
  Status RecomputeAllSegments();

  // Flows the given |segment| for the duration of |manœuvre| using its
  // intrinsic acceleration.
  Status BurnSegment(NavigationManœuvre const& manœuvre,
                     not_null<DiscreteTrajectory<Barycentric>*> segment);

  // Flows the given |segment| until |desired_final_time| with no intrinsic
  // acceleration.
  Status CoastSegment(Instant const& desired_final_time,
                      not_null<DiscreteTrajectory<Barycentric>*> segment);

  // Computes new trajectories and appends them to |segments_|.  This updates
  // the last coast of |segments_| and then appends one coast and one burn for
  // each manœuvre in |manœuvres|.  If one of the integration returns an error,
  // returns that error.  In this case the trajectories that follow the one in
  // error are of length 0 and are anomalous.
  // TODO(phl): The argument should really be an std::span, but then Apple has
  // invented the Macintosh.
  Status ComputeSegments(std::vector<NavigationManœuvre>::iterator begin,
                         std::vector<NavigationManœuvre>::iterator end);

  // Adds a trajectory to |segments_|, forked at the end of the last one.  If
  // there are already anomalous trajectories, the newly created trajectory is
  // anomalous too.
  void AddLastSegment();

  // Forgets the last trajectory after its fork.  If that trajectory was the
  // only anomalous one, there are no anomalous trajectories after this call.
  void ResetLastSegment();

  // Deletes the last trajectory and removes it from |segments_|.  If there are
  // anomalous trajectories, their number is decremented and may become 0.
  void PopLastSegment();

  Instant start_of_last_coast() const;

  // In the following functions, |index| refers to the index of a manœuvre.
  Instant start_of_next_burn(int index) const;
  Instant start_of_previous_coast(int index) const;

  Mass const initial_mass_;
  Instant initial_time_;
  DegreesOfFreedom<Barycentric> initial_degrees_of_freedom_;
  Instant desired_final_time_;
  // The root of the flight plan.  Contains a single point, not part of
  // |segments_|.  Owns all the |segments_|.
  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> root_;

  // Never empty; Starts and ends with a coast; coasts and burns alternate.
  // Each trajectory is a fork of the previous one.
  std::vector<not_null<DiscreteTrajectory<Barycentric>*>> segments_;
  // The last |anomalous_segments_| of |segments_| are anomalous, i.e. they
  // either end prematurely or follow an anomalous trajectory; in the latter
  // case they are empty.
  int anomalous_segments_ = 0;

  std::vector<NavigationManœuvre> manœuvres_;
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters_;
  Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
      generalized_adaptive_step_parameters_;
};

}  // namespace internal_flight_plan

using internal_flight_plan::FlightPlan;

}  // namespace ksp_plugin
}  // namespace principia
