
#pragma once

#include <list>
#include <memory>
#include <vector>

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

using base::IteratorOn;
using base::not_null;
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
  Vessel(not_null<Celestial const*> parent,
         not_null<Ephemeris<Barycentric>*> ephemeris,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             prediction_adaptive_step_parameters);

  Vessel(Vessel const&) = delete;
  Vessel(Vessel&&) = delete;
  Vessel& operator=(Vessel const&) = delete;
  Vessel& operator=(Vessel&&) = delete;

  virtual ~Vessel() = default;

  // Returns the body for this vessel.
  virtual not_null<MasslessBody const*> body() const;

  virtual not_null<Celestial const*> parent() const;
  virtual void set_parent(not_null<Celestial const*> parent);

  // Adds the given part to this vessel.
  virtual void AddPart(not_null<std::unique_ptr<Part>> part);
  // Removes and returns the part with the given ID.  This may empty |parts_|,
  // as happens when a vessel ceases to exist while loaded.  Note that in that
  // case |FreeParts| must not be called.
  virtual not_null<std::unique_ptr<Part>> ExtractPart(PartId id);
  // Prevents the part with the given ID from being removed in the next call to
  // |FreeParts|.
  virtual void KeepPart(PartId id);
  // Removes any part for which |AddPart| or |KeepPart| has not been called
  // since the last call to |FreePart|.  Checks that there are still parts left
  // after the removals; thus a call to |AddParts| must occur before |FreeParts|
  // is first called.
  virtual void FreeParts();

  // Returns the part with the given ID.  Such a part must have been added using
  // |AddPart|.
  virtual not_null<Part*> part(PartId id) const;

  // Calls |action| on one (non-dummy) part.
  virtual void ForSomePart(std::function<void(Part&)> action) const;
  // Calls |action| on all non-dummy parts.
  virtual void ForAllParts(std::function<void(Part&)> action) const;

  virtual DiscreteTrajectory<Barycentric> const& prediction() const;

  virtual void Vessel::set_prediction_adaptive_step_parameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          prediction_adaptive_step_parameters);
  virtual Ephemeris<Barycentric>::AdaptiveStepParameters const&
  prediction_adaptive_step_parameters() const;

  // Requires |has_flight_plan()|.
  virtual FlightPlan& flight_plan() const;
  virtual bool has_flight_plan() const;

  virtual void AdvanceTime(Instant const& time);

  // Forgets the trajectories and flight plan before |time|.  This may delete
  // the flight plan.
  virtual void ForgetBefore(Instant const& time);

  // Creates a |flight_plan_| at the end of history using the given parameters.
  // Deletes any pre-existing predictions.
  virtual void CreateFlightPlan(
      Instant const& final_time,
      Mass const& initial_mass,
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          flight_plan_adaptive_step_parameters);

  // Deletes the |flight_plan_|.  Performs no action unless |has_flight_plan()|.
  virtual void DeleteFlightPlan();

  virtual void UpdatePrediction(Instant const& last_time);

  virtual DiscreteTrajectory<Barycentric> const& psychohistory() const;
  virtual bool psychohistory_is_authoritative() const;

  // The vessel must satisfy |is_initialized()|.
  virtual void WriteToMessage(not_null<serialization::Vessel*> message) const;
  static not_null<std::unique_ptr<Vessel>> ReadFromMessage(
      serialization::Vessel const& message,
      not_null<Ephemeris<Barycentric>*> ephemeris,
      not_null<Celestial const*> parent);

 protected:
  // For mocking.
  Vessel();

 private:
  void AppendToPsychohistory(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom,
      bool authoritative);

  void FlowPrediction(Instant const& time);

  // Returns the last authoritative point of the psychohistory.
  DiscreteTrajectory<Barycentric>::Iterator last_authoritative() const;

  MasslessBody const body_;
  Ephemeris<Barycentric>::AdaptiveStepParameters
      prediction_adaptive_step_parameters_;
  // The parent body for the 2-body approximation. Not owning.
  not_null<Celestial const*> parent_;
  not_null<Ephemeris<Barycentric>*> const ephemeris_;

  std::unique_ptr<Part> dummy_part_;
  std::map<PartId, not_null<std::unique_ptr<Part>>> parts_;
  std::set<not_null<Part const*>> kept_parts_;

  // The new implementation of history, also encompasses the prolongation.
  DiscreteTrajectory<Barycentric> psychohistory_;
  bool psychohistory_is_authoritative_;

  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> prediction_;

  std::unique_ptr<FlightPlan> flight_plan_;
};

// Factories for use by the clients and the compatibility code.
Ephemeris<Barycentric>::FixedStepParameters DefaultHistoryParameters();
Ephemeris<Barycentric>::AdaptiveStepParameters DefaultProlongationParameters();
Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPredictionParameters();

}  // namespace internal_vessel

using internal_vessel::DefaultHistoryParameters;
using internal_vessel::DefaultPredictionParameters;
using internal_vessel::DefaultProlongationParameters;
using internal_vessel::Vessel;

}  // namespace ksp_plugin
}  // namespace principia
