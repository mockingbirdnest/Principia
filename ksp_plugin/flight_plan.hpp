#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "base/jthread.hpp"
#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/orbit_analyser.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan {
namespace internal {

using namespace principia::base::_jthread;
using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_orbit_analyser;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_discrete_trajectory_segment_iterator;
using namespace principia::physics::_ephemeris;
using namespace principia::quantities::_quantities;

// A chain of trajectories obtained by executing the corresponding
// `NavigationManœuvre`s.
class FlightPlan {
 public:
  // Creates a `FlightPlan` with no burns starting at `initial_time` with
  // `initial_degrees_of_freedom` and with the given `initial_mass`.  The
  // trajectories are computed using the given parameters by the given
  // `ephemeris`.  The flight plan contains a single coast which, if possible
  // ends at `desired_final_time`.
  FlightPlan(Mass const& initial_mass,
             Instant const& initial_time,
             DegreesOfFreedom<Barycentric> initial_degrees_of_freedom,
             Instant const& desired_final_time,
             not_null<Ephemeris<Barycentric>*> ephemeris,
             Ephemeris<Barycentric>::AdaptiveStepParameters
                 adaptive_step_parameters,
             Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
                 generalized_adaptive_step_parameters);

  explicit FlightPlan(FlightPlan const& other);

  virtual ~FlightPlan() = default;

  // Construction parameters.
  virtual Instant initial_time() const;
  virtual Instant desired_final_time() const;

  // End time of the last coast.  If this is less than `desired_final_time()`,
  // there is at least an anomalous manœuvre.
  virtual Instant actual_final_time() const;

  // The number of manœuvres in the flight plan.
  virtual int number_of_manœuvres() const;

  // The number of manœuvres at the end of flight plan that are anomalous, i.e.,
  // follow an anomalous segment.  These are manœuvres whose Frenet trihedron
  // cannot be drawn.  The functions that change manœuvres may change the number
  // of anomalous manœuvres.
  virtual int number_of_anomalous_manœuvres() const;

  // Returns the status associated with the first anomalous segment.  Note that
  // this may be non-OK even if `number_of_anomalous_manœuvres()` is 0, in the
  // case where only the final coast is anomalous.  The functions that change
  // manœuvres as well as the function that "avoid deadlines" may change this
  // status.
  virtual absl::Status const& anomalous_status() const;

  // Returns the specified manœuvre.  `index` must be in
  // [0, number_of_manœuvres()[.
  virtual NavigationManœuvre const& GetManœuvre(int index) const;

  // Inserts a manœuvre at the given `index` using the specified `burn`. `index`
  // must be in [0, number_of_manœuvres()].  Returns an error and has no effect
  // if the given `burn` cannot fit between the preceding burn and the following
  // one or the end of the flight plan.  Otherwise, updates the flight plan and
  // returns the integration status.
  virtual absl::Status Insert(NavigationManœuvre::Burn const& burn, int index);

  // Removes the manœuvre with the given `index`, which must be in
  // [0, number_of_manœuvres()[.
  virtual absl::Status Remove(int index);

  // Replaces a manœuvre with one using the specified `burn`.  `index` must be
  // in [0, number_of_manœuvres()[.  Returns an error and has no effect if the
  // given `burn` cannot fit between the preceding and following burns.
  // Otherwise, updates the flight plan and returns the integration status.
  virtual absl::Status Replace(NavigationManœuvre::Burn const& burn, int index);

  // Updates the desired final time of the flight plan.  Returns an error and
  // has no effect if `desired_final_time` is before the beginning of the last
  // coast.
  virtual absl::Status SetDesiredFinalTime(Instant const& desired_final_time);

  // Sets the parameters used to compute the trajectories and recomputes them
  // all.  Returns the integration status.
  virtual absl::Status SetAdaptiveStepParameters(
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

  // `index` must be in [0, number_of_segments()[.
  virtual DiscreteTrajectorySegmentIterator<Barycentric>
  GetSegment(int index) const;
  virtual DiscreteTrajectory<Barycentric> const& GetAllSegments() const;

  // Same as above, but if the flight plan is anomalous because of a deadline,
  // tries to recompute it in case the ephemeris is long enough.  This can still
  // run into a deadline.
  virtual DiscreteTrajectorySegmentIterator<Barycentric>
  GetSegmentAvoidingDeadlines(int index);
  virtual DiscreteTrajectory<Barycentric> const&
  GetAllSegmentsAvoidingDeadlines();

  // Orbit analysis is enabled at construction, and may be enabled/disabled
  // dynamically.
  void EnableAnalysis(bool enabled);

  // `coast_index` must be in [0, number_of_manœuvres()].
  virtual OrbitAnalyser::Analysis* analysis(int coast_index);
  double progress_of_analysis(int coast_index) const;

  void WriteToMessage(not_null<serialization::FlightPlan*> message) const;

  // This may return a null pointer if the flight plan contained in the
  // `message` is anomalous.
  static std::unique_ptr<FlightPlan> ReadFromMessage(
      serialization::FlightPlan const& message,
      not_null<Ephemeris<Barycentric>*> ephemeris);

  static constexpr std::int64_t max_ephemeris_steps_per_frame = 1000;

  static constexpr absl::StatusCode bad_desired_final_time =
      absl::StatusCode::kOutOfRange;
  static constexpr absl::StatusCode does_not_fit =
      absl::StatusCode::kOutOfRange;
  static constexpr absl::StatusCode singular =
      absl::StatusCode::kInvalidArgument;

 protected:
  // For mocking.
  FlightPlan();

 private:
  // Clears and recomputes all trajectories in `segments_`.
  absl::Status RecomputeAllSegments();

  // If the flight plan is anomalous because of an integration deadline, try to
  // recompute it from the first anomalous segment.  This might work better if
  // the ephemeris has been prolonged enough.
  absl::Status RecomputeSegmentsAvoidingDeadlineIfNeeded();

  // Flows the given `segment` for the duration of `manœuvre` using its
  // intrinsic acceleration.
  absl::Status BurnSegment(
      NavigationManœuvre const& manœuvre,
      DiscreteTrajectorySegmentIterator<Barycentric> segment,
      std::int64_t max_ephemeris_steps);

  // Flows the given `segment` until `desired_final_time` with no intrinsic
  // acceleration.
  absl::Status CoastSegment(
      Instant const& desired_final_time,
      DiscreteTrajectorySegmentIterator<Barycentric> segment,
      std::int64_t max_ephemeris_steps);

  // Computes new trajectories and appends them to `segments_`.  This updates
  // the last coast of `segments_` and then appends one coast and one burn for
  // each manœuvre in `manœuvres`.  If one of the integration returns an error,
  // returns that error.  In this case the trajectories that follow the one in
  // error are of length 0 and are anomalous.
  // TODO(phl): The argument should really be an std::span, but then Apple has
  // invented the Macintosh.
  absl::Status ComputeSegments(std::vector<NavigationManœuvre>::iterator begin,
                               std::vector<NavigationManœuvre>::iterator end,
                               std::int64_t max_ephemeris_steps);

  // Adds a trajectory to `segments_`, forked at the end of the last one.  If
  // there are already anomalous trajectories, the newly created trajectory is
  // anomalous too.
  void AddLastSegment();

  // Forgets the last trajectory after its fork.  If that trajectory was the
  // only anomalous one, there are no anomalous trajectories after this call.
  void ResetLastSegment();

  // Deletes the last `count` trajectories and removes them from `segments_`. If
  // there are anomalous trajectories, their number is decremented and may
  // become 0.  `count` must be even (possibly 0) to maintain the invariant that
  // the number of segments is odd.
  void PopLastSegments(std::int64_t count);

  // Pops the burn of the manœuvre with the given index and all following
  // segments, then resets the last segment (which is the coast preceding
  // `manœuvres_[index]`).
  void PopSegmentsAffectedByManœuvre(int index);

  // Reconstructs each manœuvre after `manœuvres_[index]` (starting with
  // `manœuvres_[index + 1]`), keeping the same burns but recomputing the
  // initial masses from `manœuvres_[index].final_mass()`.
  void UpdateInitialMassOfManœuvresAfter(int index);

  // Starts a thread to prolong the ephemeris if needed.
  void MakeProlongator(Instant const& prolongation_time);

  Instant start_of_last_coast() const;

  // In the following functions, `index` refers to the index of a manœuvre.
  Instant start_of_burn(int index) const;
  Instant start_of_next_burn(int index) const;
  Instant start_of_previous_coast(int index) const;

  Mass const initial_mass_;
  Instant const initial_time_;
  DegreesOfFreedom<Barycentric> const initial_degrees_of_freedom_;
  not_null<Ephemeris<Barycentric>*> const ephemeris_;

  Instant desired_final_time_;

  // The trajectory of the part, composed of any number of segments,
  // alternatively coasts and burns.
  DiscreteTrajectory<Barycentric> trajectory_;

  // Never empty; Starts and ends with a coast; coasts and burns alternate.
  std::vector<DiscreteTrajectorySegmentIterator<Barycentric>> segments_;
  // The last `anomalous_segments_` of `segments_` are anomalous, i.e., they
  // either end prematurely or follow an anomalous segment; in the latter case
  // they are empty.
  int anomalous_segments_ = 0;
  // The status of the first anomalous segment.  Set and used exclusively by
  // `ComputeSegments`.
  absl::Status anomalous_status_;

  std::vector<NavigationManœuvre> manœuvres_;

  // The coast analysers always exist for each coast, but they may be idle if
  // analysis is disabled.
  std::vector<not_null<std::unique_ptr<OrbitAnalyser>>> coast_analysers_;
  bool analysis_is_enabled_ = true;

  // These members are only accessed by the main thread.
  jthread prolongator_;
  Instant last_prolongation_time_ = InfinitePast;

  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters_;
  Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
      generalized_adaptive_step_parameters_;
};

}  // namespace internal

using internal::FlightPlan;

}  // namespace _flight_plan
}  // namespace ksp_plugin
}  // namespace principia
