
#pragma once

#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/synchronization/mutex.h"
#include "base/jthread.hpp"
#include "base/recurring_thread.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/orbit_analyser.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "physics/checkpointer.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {

class VesselTest;

namespace internal_vessel {

using base::not_null;
using base::RecurringThread;
using geometry::InfinitePast;
using geometry::Instant;
using geometry::Vector;
using physics::Checkpointer;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::DiscreteTrajectorySegment;
using physics::DiscreteTrajectorySegmentIterator;
using physics::Ephemeris;
using physics::MasslessBody;
using physics::RotatingBody;
using quantities::Force;
using quantities::GravitationalParameter;
using quantities::Mass;
using quantities::Time;

// Represents a KSP |Vessel|.
class Vessel {
 public:
  using Manœuvres = std::vector<
      not_null<std::unique_ptr<Manœuvre<Barycentric, Navigation> const>>>;

  // Constructs a vessel whose parent is initially |*parent|.
  Vessel(GUID guid,
         std::string name,
         not_null<Celestial const*> parent,
         not_null<Ephemeris<Barycentric>*> ephemeris,
         Ephemeris<Barycentric>::AdaptiveStepParameters
             prediction_adaptive_step_parameters,
         DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters const&
             downsampling_parameters);

  Vessel(Vessel const&) = delete;
  Vessel(Vessel&&) = delete;
  Vessel& operator=(Vessel const&) = delete;
  Vessel& operator=(Vessel&&) = delete;

  virtual ~Vessel();

  // Returns the GUID passed at construction.
  virtual GUID const& guid() const;

  // Returns the name.
  virtual std::string const& name() const;
  // Changes the name.
  virtual void set_name(std::string const& new_name);

  // Returns the body for this vessel.
  virtual not_null<MasslessBody const*> body() const;

  virtual not_null<Celestial const*> parent() const;
  virtual void set_parent(not_null<Celestial const*> parent);

  // Adds the given part to this vessel.  Note that this does not add the part
  // to the set of kept parts, and that unless |KeepPart| is called, the part
  // will be removed by the next call to |FreeParts|.
  virtual void AddPart(not_null<std::unique_ptr<Part>> part);
  // Removes and returns the part with the given ID.  This may empty |parts_|,
  // as happens when a vessel ceases to exist while loaded.  Note that in that
  // case |FreeParts| must not be called.
  virtual not_null<std::unique_ptr<Part>> ExtractPart(PartId id);
  // Prevents the part with the given ID from being removed in the next call to
  // |FreeParts|.
  virtual void KeepPart(PartId id);
  // Whether |KeepPart| was called with this |id| since the last call to
  // |FreeParts|.
  virtual bool WillKeepPart(PartId id) const;
  // Removes any part for which |KeepPart| has not been called since the last
  // call to |FreeParts|.  Checks that there are still parts left after the
  // removals; thus a call to |AddPart| must occur before |FreeParts| is first
  // called.
  virtual void FreeParts();

  // Clears the forces and torques on all parts.
  virtual void ClearAllIntrinsicForcesAndTorques();

  // Detects a change in the collapsibility of the vessel and creates a new
  // trajectory segment if needed.  Must be called after the pile-ups have been
  // collected.
  virtual void DetectCollapsibilityChange() EXCLUDES(lock_);

  // If the history is empty, appends a single point to it, computed as the
  // barycentre of all parts.  |parts_| must not be empty.  After this call,
  // |history_| is never empty again and the psychohistory is usable.  Must be
  // called (at least once) after the creation of the vessel.
  virtual void CreateHistoryIfNeeded(Instant const& t) EXCLUDES(lock_);

  // Disables downsampling for the history of this vessel.  This is useful when
  // the vessel collided with a celestial, as downsampling might run into
  // trouble.
  virtual void DisableDownsampling() EXCLUDES(lock_);

  // Returns the part with the given ID.  Such a part must have been added using
  // |AddPart|.
  virtual not_null<Part*> part(PartId id) const;

  // Calls |action| on one part.
  virtual void ForSomePart(std::function<void(Part&)> action) const;
  // Calls |action| on all parts.
  virtual void ForAllParts(std::function<void(Part&)> action) const;

  virtual DiscreteTrajectory<Barycentric> const& trajectory() const;
  virtual DiscreteTrajectorySegmentIterator<Barycentric> history() const;
  virtual DiscreteTrajectorySegmentIterator<Barycentric> psychohistory() const;
  virtual DiscreteTrajectorySegmentIterator<Barycentric> prediction() const;

  virtual void set_prediction_adaptive_step_parameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          prediction_adaptive_step_parameters);
  virtual Ephemeris<Barycentric>::AdaptiveStepParameters const&
  prediction_adaptive_step_parameters() const;

  // Requires |has_flight_plan()|.
  virtual FlightPlan& flight_plan() const;
  virtual bool has_flight_plan() const;

  // Extends the history and psychohistory of this vessel by computing the
  // centre of mass of its parts at every point in their history and
  // psychohistory.  Clears the parts' history and psychohistory.
  virtual void AdvanceTime();

  // Asks the reanimator thread to asynchronously reconstruct the past so that
  // the |t_min()| of the vessel ultimately ends up at or before
  // |desired_t_min|.
  void RequestReanimation(Instant const& desired_t_min);

  // Blocks until the |t_min()| of the vessel is at or before |desired_t_min|.
  void WaitForReanimation(Instant const& desired_t_min);

  // Creates a |flight_plan_| at the end of history using the given parameters.
  virtual void CreateFlightPlan(
      Instant const& final_time,
      Mass const& initial_mass,
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          flight_plan_adaptive_step_parameters,
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
          flight_plan_generalized_adaptive_step_parameters) EXCLUDES(lock_);

  // Deletes the |flight_plan_|.  Performs no action unless |has_flight_plan()|.
  virtual void DeleteFlightPlan();

  // Requires |has_flight_plan()|.
  // If |history_->back().time| lies within a planned manœuvre, UNAVAILABLE is
  // returned.
  // Otherwise, deletes |flight_plan_| and recreates it from the current
  // |history_| and the given |initial_mass|, re-adding any future manœuvres.
  // Past manœuvres are discarded, under the assumption that they have been
  // performed.
  // If |history_->back().time| is greater than the current desired final time,
  // the flight plan length is kept; otherwise, the desired final time is kept.
  absl::Status RebaseFlightPlan(Mass const& initial_mass) EXCLUDES(lock_);

  // Tries to replace the current prediction with a more recently computed one.
  // No guarantees that this happens.  No guarantees regarding the end time of
  // the prediction when this call returns.
  virtual void RefreshPrediction();

  // Same as above, but when this call returns the prediction is guaranteed to
  // have a last time at or before |time|.
  virtual void RefreshPrediction(Instant const& time);

  // Stop the asynchronous prognosticator as soon as convenient.
  void StopPrognosticator();

  // Stops any analyser running for a different mission duration and triggers a
  // new analysis.
  void RequestOrbitAnalysis(Time const& mission_duration);

  // Stops the analyser.
  void ClearOrbitAnalyser();

  // Returns a number between 0 and 1 indicating how far we are within the
  // current analysis.
  double progress_of_orbit_analysis() const;

  // Prepares the last completed analysis so that will be returned by
  // |orbit_analysis|.
  // TODO(phl): This API is weird.  Why does the caller need a 2-step dance?
  void RefreshOrbitAnalysis();

  // Returns the latest completed analysis, if there is one.
  OrbitAnalyser::Analysis* orbit_analysis();

  // Returns "vessel_name (GUID)".
  std::string ShortDebugString() const;

  // The vessel must satisfy |is_initialized()|.
  virtual void WriteToMessage(not_null<serialization::Vessel*> message,
                              PileUp::SerializationIndexForPileUp const&
                                  serialization_index_for_pile_up) const
      EXCLUDES(lock_);
  static not_null<std::unique_ptr<Vessel>> ReadFromMessage(
      serialization::Vessel const& message,
      not_null<Celestial const*> parent,
      not_null<Ephemeris<Barycentric>*> ephemeris,
      std::function<void(PartId)> const& deletion_callback);
  void FillContainingPileUpsFromMessage(
      serialization::Vessel const& message,
      PileUp::PileUpForSerializationIndex const&
          pile_up_for_serialization_index);

  static void MakeAsynchronous();
  static void MakeSynchronous();

 protected:
  // For mocking.
  Vessel();

 private:
  struct PrognosticatorParameters {
    Instant first_time;
    DegreesOfFreedom<Barycentric> first_degrees_of_freedom;
    Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters;
  };
  friend bool operator!=(PrognosticatorParameters const& left,
                         PrognosticatorParameters const& right);

  using TrajectoryIterator =
      DiscreteTrajectory<Barycentric>::iterator (Part::*)();

  // Return functions that can be passed to a |Checkpointer| to write this
  // vessel to a checkpoint or read it back.
  Checkpointer<serialization::Vessel>::Writer
  MakeCheckpointerWriter();
  Checkpointer<serialization::Vessel>::Reader
  MakeCheckpointerReader();

  absl::Status Reanimate(Instant const desired_t_min) EXCLUDES(lock_);

  // |t_initial| is the time of the checkpoint, which is the end of the non-
  // collapsible segment.
  absl::Status ReanimateOneCheckpoint(
      serialization::Vessel::Checkpoint const& message,
      Instant const& t_initial) EXCLUDES(lock_);

  // Runs the integrator to compute the |prognostication_| based on the given
  // parameters.
  absl::StatusOr<DiscreteTrajectory<Barycentric>>
  FlowPrognostication(PrognosticatorParameters prognosticator_parameters);

  // Appends to |trajectory_| the centre of mass of the trajectories of the
  // parts denoted by |part_trajectory_begin| and |part_trajectory_end|.  Only
  // the points that are strictly after the start of the |segment| are used.
  void AppendToVesselTrajectory(
      TrajectoryIterator part_trajectory_begin,
      TrajectoryIterator part_trajectory_end,
      DiscreteTrajectorySegment<Barycentric> const& segment);

  // Attaches the given |trajectory| to the end of the |psychohistory_| to
  // become the new |prediction_|.  If |prediction_| is not null, it is deleted.
  void AttachPrediction(DiscreteTrajectory<Barycentric>&& trajectory);

  // A vessel is collapsible if it is alone in its pile-up and is in inertial
  // motion.
  bool IsCollapsible() const;

  GUID const guid_;
  std::string name_;

  MasslessBody const body_;
  Ephemeris<Barycentric>::AdaptiveStepParameters
      prediction_adaptive_step_parameters_;
  // The parent body for the 2-body approximation.
  not_null<Celestial const*> parent_;
  not_null<Ephemeris<Barycentric>*> const ephemeris_;
  // TODO(phl): Serialize this field.
  DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters const
      downsampling_parameters_;

  // The techniques and terminology follow [Lov22].
  RecurringThread<Instant> reanimator_;

  // TODO(phl): Verify locking.
  mutable absl::Mutex lock_;

  // When reading a pre-Zermelo save, the existing history must be
  // non-collapsible as we don't know anything about it.
  bool is_collapsible_ = false;

  std::map<PartId, not_null<std::unique_ptr<Part>>> parts_;
  std::set<PartId> kept_parts_;

  // The vessel trajectory is made of the history (always present) and (most of
  // the time) the psychohistory and prediction.  The prediction is periodically
  // recomputed by the prognosticator.
  DiscreteTrajectory<Barycentric> trajectory_;

  not_null<std::unique_ptr<Checkpointer<serialization::Vessel>>> checkpointer_;

  // This member must only be accessed by the |reanimator_| thread, or before
  // the |reanimator_| thread is started.
  Instant oldest_reanimated_checkpoint_ = InfinitePast;

  // Parameter passed to the last call to |RequestReanimation|, if any.
  std::optional<Instant> last_desired_t_min_ GUARDED_BY(lock_);

  // See the comments in pile_up.hpp for an explanation of the terminology.
  // The |history_| is empty until the first call to CreateHistoryIfNeeded.
  // It is made of a series of forks, alternatively non-collapsible and
  // collapsible.
  DiscreteTrajectorySegmentIterator<Barycentric> history_
      GUARDED_BY(lock_);

  // The last (most recent) segment of the |history_| prior to the
  // |psychohistory_|.  May be identical to |history_|.  Always identical to
  // |std::prev(psychohistory_)|.
  DiscreteTrajectorySegmentIterator<Barycentric> backstory_;

  // The |psychohistory_| is forked off the end of the |history_| and the
  // |prediction_| is forked off the end of the |psychohistory_|.
  DiscreteTrajectorySegmentIterator<Barycentric> psychohistory_;
  DiscreteTrajectorySegmentIterator<Barycentric> prediction_;

  RecurringThread<PrognosticatorParameters,
                  DiscreteTrajectory<Barycentric>> prognosticator_;

  std::unique_ptr<FlightPlan> flight_plan_;

  std::optional<OrbitAnalyser> orbit_analyser_;

  static std::atomic_bool synchronous_;

  friend class ksp_plugin::VesselTest;
};

}  // namespace internal_vessel

using internal_vessel::Vessel;

}  // namespace ksp_plugin
}  // namespace principia
