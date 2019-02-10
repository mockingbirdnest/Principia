
#pragma once

#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "base/status.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using base::not_null;
using base::Status;
using geometry::Instant;
using geometry::Vector;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MasslessBody;
using quantities::Force;
using quantities::GravitationalParameter;
using quantities::Mass;

// Represents a KSP |Vessel|.
class Vessel {
 public:
  using Manœuvres = std::vector<
      not_null<std::unique_ptr<Manœuvre<Barycentric, Navigation> const>>>;

  // Constructs a vessel whose parent is initially |*parent|.  No transfer of
  // ownership.
  Vessel(GUID const& guid,
         std::string const& name,
         not_null<Celestial const*> parent,
         not_null<Ephemeris<Barycentric>*> ephemeris,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             prediction_adaptive_step_parameters);

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

  virtual void ClearAllIntrinsicForces();

  // If the history is empty, appends a single point to it, computed as the
  // barycentre of all parts.  |parts_| must not be empty.  After this call,
  // |history_| is never empty again and the psychohistory is usable.
  virtual void PrepareHistory(Instant const& t);

  // Disables downsampling for the history of this vessel.  This is useful when
  // the vessel collided with a celestial, as downsampling might run into
  // trouble.
  virtual void DisableDownsampling();

  // Returns the part with the given ID.  Such a part must have been added using
  // |AddPart|.
  virtual not_null<Part*> part(PartId id) const;

  // Calls |action| on one part.
  virtual void ForSomePart(std::function<void(Part&)> action) const;
  // Calls |action| on all parts.
  virtual void ForAllParts(std::function<void(Part&)> action) const;

  virtual DiscreteTrajectory<Barycentric> const& psychohistory() const;
  virtual DiscreteTrajectory<Barycentric> const& prediction() const;

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

  // Forgets the trajectories and flight plan before |time|.  This may delete
  // the flight plan.
  virtual void ForgetBefore(Instant const& time);

  // Creates a |flight_plan_| at the end of history using the given parameters.
  virtual void CreateFlightPlan(
      Instant const& final_time,
      Mass const& initial_mass,
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          flight_plan_adaptive_step_parameters,
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
          flight_plan_generalized_adaptive_step_parameters);

  // Deletes the |flight_plan_|.  Performs no action unless |has_flight_plan()|.
  virtual void DeleteFlightPlan();

  // Tries to extend the prediction by extending the ephemeris by at most
  // |max_ephemeris_steps_per_frame|.  No guarantees regarding the end time of
  // the prediction when this call returns.
  virtual void FlowPrediction();

  // Extends the prediction (and the ephemeris) up to and including |time|.  May
  // not be able to do so next to a singularity, in which case an error is
  // returned.
  virtual Status FlowPrediction(Instant const& time);

  // The vessel must satisfy |is_initialized()|.
  virtual void WriteToMessage(not_null<serialization::Vessel*> message,
                              PileUp::SerializationIndexForPileUp const&
                                  serialization_index_for_pile_up) const;
  static not_null<std::unique_ptr<Vessel>> ReadFromMessage(
      serialization::Vessel const& message,
      not_null<Celestial const*> parent,
      not_null<Ephemeris<Barycentric>*> ephemeris,
      std::function<void(PartId)> const& deletion_callback);
  void FillContainingPileUpsFromMessage(
      serialization::Vessel const& message,
      PileUp::PileUpForSerializationIndex const&
          pile_up_for_serialization_index);

  // Returns "vessel_name (GUID)".
  std::string ShortDebugString() const;

 protected:
  // For mocking.
  Vessel();

 private:
  struct PrognosticatorParameters {
    Instant first_time;
    DegreesOfFreedom<Barycentric> first_degrees_of_freedom;
    std::optional<Instant> last_time;
    Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters;
    bool shutdown = false;
  };

  using TrajectoryIterator =
      DiscreteTrajectory<Barycentric>::Iterator (Part::*)();

  // Run by the |prognosticator_| thread to periodically recompute the
  // prognostication.
  void RepeatedlyFlowPrognostication();

  // Appends to |trajectory| the centre of mass of the trajectories of the parts
  // denoted by |part_trajectory_begin| and |part_trajectory_end|.
  void AppendToVesselTrajectory(TrajectoryIterator part_trajectory_begin,
                                TrajectoryIterator part_trajectory_end,
                                DiscreteTrajectory<Barycentric>& trajectory);

  // Attaches the given |trajectory| to the end of the |psychohistory_| to
  // become the new |prediction_|.
  void AttachPrediction(
      not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> trajectory);

  GUID const guid_;
  std::string name_;

  MasslessBody const body_;
  Ephemeris<Barycentric>::AdaptiveStepParameters
      prediction_adaptive_step_parameters_;
  // The parent body for the 2-body approximation.
  not_null<Celestial const*> parent_;
  not_null<Ephemeris<Barycentric>*> const ephemeris_;

  std::map<PartId, not_null<std::unique_ptr<Part>>> parts_;
  std::set<PartId> kept_parts_;

  mutable absl::Mutex prognosticator_lock_;
  std::optional<PrognosticatorParameters> prognosticator_parameters_
      GUARDED_BY(prognosticator_lock_);
  Status prognosticator_status_ GUARDED_BY(prognosticator_lock_);
  std::thread prognosticator_;

  // See the comments in pile_up.hpp for an explanation of the terminology.
  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> history_;
  DiscreteTrajectory<Barycentric>* psychohistory_ = nullptr;

  // The |prediction_| is forked off the end of the |psychohistory_|.
  DiscreteTrajectory<Barycentric>* prediction_ = nullptr;

  // The |prognostication_| is a root trajectory that's computed asynchronously
  // and may or may not be used as a prediction;
  std::unique_ptr<DiscreteTrajectory<Barycentric>> prognostication_
      GUARDED_BY(prognosticator_lock_);

  std::unique_ptr<FlightPlan> flight_plan_;
};

}  // namespace internal_vessel

using internal_vessel::Vessel;

}  // namespace ksp_plugin
}  // namespace principia
