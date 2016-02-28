
#pragma once

#include <memory>
#include <vector>

#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/part.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MasslessBody;
using quantities::GravitationalParameter;

namespace ksp_plugin {

// Represents a KSP |Vessel|.
class Vessel {
 public:
  using Manœuvres =
      std::vector<
          not_null<std::unique_ptr<Manœuvre<Barycentric, Navigation> const>>>;

  Vessel(Vessel const&) = delete;
  Vessel(Vessel&&) = delete;
  Vessel& operator=(Vessel const&) = delete;
  Vessel& operator=(Vessel&&) = delete;
  ~Vessel() = default;

  // Constructs a vessel whose parent is initially |*parent|.  No transfer of
  // ownership.
  explicit Vessel(not_null<Celestial const*> const parent);

  // Returns the body for this vessel.
  virtual not_null<MasslessBody const*> body() const;

  // True if, and only if, |prolongation_| is not null, i.e., if either
  // |CreateProlongation| or |CreateHistoryAndForkProlongation| was called at
  // some point.
  virtual bool is_initialized() const;

  virtual not_null<Celestial const*> parent() const;
  virtual void set_parent(not_null<Celestial const*> const parent);

  virtual DiscreteTrajectory<Barycentric> const& history() const;

  // Both accessors require |is_initialized()|.
  virtual DiscreteTrajectory<Barycentric> const& prolongation() const;
  virtual not_null<DiscreteTrajectory<Barycentric>*> mutable_prolongation();

  // Requires |is_initialized()|.
  virtual not_null<FlightPlan*> flight_plan() const;
  virtual bool has_flight_plan() const;

  // Requires |has_prediction()|.
  virtual DiscreteTrajectory<Barycentric> const& prediction() const;
  virtual bool has_prediction() const;

  // Creates a |history_| for this vessel and appends a point with the
  // given |time| and |degrees_of_freedom|, then forks a |prolongation_| at
  // |time|.  Nulls |owned_prolongation_|.  The vessel must not satisfy
  // |is_initialized()|.  The vessel |is_initialized()| after the call.
  virtual void CreateHistoryAndForkProlongation(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Appends a point to the history.
  virtual void AppendToHistory(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  //TODO(phl):comment
  virtual void AdvanceTime(Instant const& time);

  // Forgets the history before |time|.  |time| must not be after
  // |ForgettableTime|.
  virtual void ForgetBefore(Instant const& time);

  // Returns the largest time for which |ForgetBefore| may be called.
  virtual Instant ForgettableTime() const;

  // Creates a |flight_plan_| at the end of history using the given parameters.
  // Deletes any pre-existing predictions.
  virtual void CreateFlightPlan(
      Instant const& final_time,
      Mass const& initial_mass,
      not_null<Ephemeris<Barycentric>*> ephemeris,
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          adaptive_parameters);

  // Deletes the |flight_plan_|.  Performs no action unless |has_flight_plan()|.
  virtual void DeleteFlightPlan();

  virtual void UpdatePrediction(
      not_null<Ephemeris<Barycentric>*> ephemeris,
      Instant const& last_time,
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          adaptive_parameters);

  // Deletes the |prediction_|.  Performs no action unless |has_prediction()|.
  virtual void DeletePrediction();

  // The vessel must satisfy |is_initialized()|.
  virtual void WriteToMessage(
      not_null<serialization::Vessel*> const message) const;
  static not_null<std::unique_ptr<Vessel>> ReadFromMessage(
      serialization::Vessel const& message,
      not_null<Ephemeris<Barycentric>*> const ephemeris,
      not_null<Celestial const*> const parent);

 protected:
  // For mocking.
  Vessel();

 private:
  // Deletes the |prolongation_| and forks a new one at |time|.
  // The vessel must satisfy |is_initialized()|.
  void ResetProlongation(Instant const& time);

  MasslessBody const body_;
  // The parent body for the 2-body approximation. Not owning.
  not_null<Celestial const*> parent_;
  // The past and present trajectory of the body. It ends at |HistoryTime()|
  // unless |*this| was created after |HistoryTime()|, in which case it ends
  // at |current_time_|.  It is advanced with a constant time step.
  std::unique_ptr<DiscreteTrajectory<Barycentric>> history_;
  // Most of the time, this is a child trajectory of |*history_|. It is forked
  // at |history_->last_time()| and continues until |current_time_|. It is
  // computed with a non-constant timestep, which breaks symplecticity.
  // If |history_| is null, this points to |owned_prolongation_| instead.
  // Not owning.
  DiscreteTrajectory<Barycentric>* prolongation_ = nullptr;
  // Child trajectory of |history_|.
  DiscreteTrajectory<Barycentric>* prediction_ = nullptr;

  std::unique_ptr<FlightPlan> flight_plan_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/vessel_body.hpp"
