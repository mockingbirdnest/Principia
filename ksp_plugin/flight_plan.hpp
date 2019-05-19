
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
  // trajectories are computed using the given |integrator| in the given
  // |ephemeris|.
  //TODO(phl): 1 coast till desired final time.
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

  virtual Instant initial_time() const;
  virtual Instant actual_final_time() const;
  virtual Instant desired_final_time() const;

  virtual int number_of_manœuvres() const;
  //TODO(phl):comment
  virtual int number_of_anomalous_manœuvres() const;

  // |index| must be in [0, number_of_manœuvres()[.
  virtual NavigationManœuvre const& GetManœuvre(int index) const;

  // The following function returns an error and has no effect if the given
  // |burn| would start before |initial_time_| or before the end of the previous
  // burn, or end after |desired_final_time_|, or if the |burn| is singular.
  virtual Status Append(NavigationManœuvre::Burn const& burn);

  // Forgets the flight plan at least before |time|.  The actual cutoff time
  // will be in a coast trajectory and may be after |time|.  |on_empty| is run
  // if the flight plan would become empty (it is not modified before running
  // |on_empty|).
  virtual void ForgetBefore(Instant const& time,
                            std::function<void()> const& on_empty);


  // |size()| must be greater than 0.
  virtual Status RemoveLast();
  //TODO(phl):Comment.
  virtual Status Replace(NavigationManœuvre::Burn const& burn, int index);
  // |size()| must be greater than 0.
  virtual Status ReplaceLast(NavigationManœuvre::Burn const& burn);

  // Returns an error and has no effect if |desired_final_time| is before the
  // end of the last manœuvre or before |initial_time_|.
  virtual Status SetDesiredFinalTime(Instant const& desired_final_time);

  virtual Ephemeris<Barycentric>::AdaptiveStepParameters const&
  adaptive_step_parameters() const;
  virtual Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
  generalized_adaptive_step_parameters() const;

  // Sets the parameters used to compute the trajectories.  The trajectories are
  // recomputed.  Returns false (and doesn't change this object) if the
  // parameters would make it impossible to recompute the trajectories.
  virtual Status SetAdaptiveStepParameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          adaptive_step_parameters,
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
          generalized_adaptive_step_parameters);

  // Returns the number of trajectory segments in this object.
  virtual int number_of_segments() const;

  // |index| must be in [0, number_of_segments()[.  Sets the iterators to denote
  // the given trajectory segment.
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

  //TODO(phl):fix
  static constexpr Error bad_desired_final_time = Error::DEADLINE_EXCEEDED;
  static constexpr Error does_not_fit = Error::INVALID_ARGUMENT;
  static constexpr Error singular = Error::OUT_OF_RANGE;

 protected:
  // For mocking.
  FlightPlan();

 private:
  // Clears and recomputes all trajectories in |segments_|.
  Status RecomputeAllSegments();

  // Flows the last segment for the duration of |manœuvre| using its intrinsic
  // acceleration.
  //TODO(phl):comment
  Status BurnSegment(NavigationManœuvre const& manœuvre,
                     not_null<DiscreteTrajectory<Barycentric>*> segment);
  // Flows the last segment until |desired_final_time| with no intrinsic
  // acceleration.
  //TODO(phl):comment
  Status CoastSegment(Instant const& desired_final_time,
                      not_null<DiscreteTrajectory<Barycentric>*> segment);
  // TODO(phl): The first argument should really be an std::span, but then Apple
  // has invented the Macintosh.
  Status RecomputeSegments(std::vector<NavigationManœuvre>& manœuvres);

  // Adds a trajectory to |segments_|, forked at the end of the last one.
  //TODO(phl):comment
  void AddLastSegment();

  // Forgets the last segment after its fork.
  void ResetLastSegment();

  // Deletes the last segment and removes it from |segments_|.
  void PopLastSegment();

  // If the integration of a coast from the fork of |coast| until
  // |manœuvre.initial_time()| reaches the end, returns the integrated
  // trajectory.  Otherwise, returns null.
  DiscreteTrajectory<Barycentric>* CoastIfReachesManœuvreInitialTime(
      DiscreteTrajectory<Barycentric>& coast,
      NavigationManœuvre const& manœuvre);

  Instant start_of_last_coast() const;
  Instant start_of_penultimate_coast() const;

  //TODO(phl):comment
  Instant start_of_next_burn(int index) const;
  Instant start_of_previous_coast(int index) const;

  DiscreteTrajectory<Barycentric>& last_coast();
  DiscreteTrajectory<Barycentric>& penultimate_coast();

  //TODO(phl):comment
  not_null<DiscreteTrajectory<Barycentric>*> previous_coast(int index);

  Mass const initial_mass_;
  Instant initial_time_;
  DegreesOfFreedom<Barycentric> initial_degrees_of_freedom_;
  Instant desired_final_time_;
  // The root of the flight plan.  Contains a single point, not part of
  // |segments_|.  Owns all the |segments_|.
  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> root_;
  // Never empty; Starts and ends with a coasting segment; coasting and burning
  // alternate.  This simulates a stack.  Each segment is a fork of the previous
  // one.
  std::vector<not_null<DiscreteTrajectory<Barycentric>*>> segments_;
  std::vector<NavigationManœuvre> manœuvres_;
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters_;
  Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
      generalized_adaptive_step_parameters_;
  // The last |anomalous_segments_| of |segments_| are anomalous, i.e. they
  // either end prematurely or follow an anomalous segment; in the latter case
  // they are empty.
  // The contract of |Append| and |ReplaceLast| implies that
  // |anomalous_segments_| is at most 2: the penultimate coast is never
  // anomalous.
  int anomalous_segments_ = 0;
};

}  // namespace internal_flight_plan

using internal_flight_plan::FlightPlan;

}  // namespace ksp_plugin
}  // namespace principia
