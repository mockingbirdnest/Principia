
#pragma once

#include <list>
#include <memory>
#include <vector>

#include "base/disjoint_sets.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/vessel_subsets.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using base::not_null;
using base::Subset;
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
  using Manœuvres =
      std::vector<
          not_null<std::unique_ptr<Manœuvre<Barycentric, Navigation> const>>>;

  Vessel(Vessel const&) = delete;
  Vessel(Vessel&&) = delete;
  Vessel& operator=(Vessel const&) = delete;
  Vessel& operator=(Vessel&&) = delete;

  // |CHECK|s that |*this| is not piled up.
  virtual ~Vessel();

  // Constructs a vessel whose parent is initially |*parent|.  No transfer of
  // ownership.
  Vessel(not_null<Celestial const*> parent,
         not_null<Ephemeris<Barycentric>*> ephemeris,
         Ephemeris<Barycentric>::FixedStepParameters const&
             history_fixed_step_parameters,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             prolongation_adaptive_step_parameters,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             prediction_adaptive_step_parameters);

  // Returns the body for this vessel.
  virtual not_null<MasslessBody const*> body() const;

  // True if, and only if, |prolongation_| is not null, i.e., if either
  // |CreateProlongation| or |CreateHistoryAndForkProlongation| was called at
  // some point.
  virtual bool is_initialized() const;

  virtual not_null<Celestial const*> parent() const;
  virtual void set_parent(not_null<Celestial const*> parent);

  // These three functions require |is_initialized()|.
  virtual DiscreteTrajectory<Barycentric> const& history() const;
  virtual DiscreteTrajectory<Barycentric> const& prolongation() const;
  virtual DiscreteTrajectory<Barycentric> const& prediction() const;

  // Requires |has_flight_plan()|.
  virtual FlightPlan& flight_plan() const;
  virtual bool has_flight_plan() const;

  // A vessel that was in the physics since the last time its history advanced
  // (which with the time when its prolongation was reset).  For such a vessel,
  // the prolongation, not the history, is authoritative.
  virtual void set_dirty();
  virtual bool is_dirty() const;

  virtual void set_prediction_adaptive_step_parameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          prediction_adaptive_step_parameters);
  virtual Ephemeris<Barycentric>::AdaptiveStepParameters const&
      prediction_adaptive_step_parameters() const;

  // Creates a |history_| for this vessel and appends a point with the
  // given |time| and |degrees_of_freedom|, then forks a |prolongation_| at
  // |time|.  Nulls |owned_prolongation_|.  The vessel must not satisfy
  // |is_initialized()|.  The vessel |is_initialized()| after the call.
  virtual void CreateHistoryAndForkProlongation(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Advances time for a vessel not in the physics bubble.  This may clean the
  // vessel.
  virtual void AdvanceTimeNotInBubble(Instant const& time);

  // Advances time for a vessel in the physics bubble.  This dirties the vessel.
  virtual void AdvanceTimeInBubble(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

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

  virtual void AppendToPsychohistory(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom,
      bool authoritative);

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
  void AdvanceHistoryIfNeeded(Instant const& time);
  void FlowHistory(Instant const& time);
  void FlowProlongation(Instant const& time);
  void FlowPrediction(Instant const& time);

  MasslessBody const body_;
  Ephemeris<Barycentric>::FixedStepParameters const
      history_fixed_step_parameters_;
  Ephemeris<Barycentric>::AdaptiveStepParameters const
      prolongation_adaptive_step_parameters_;
  Ephemeris<Barycentric>::AdaptiveStepParameters
      prediction_adaptive_step_parameters_;
  // The parent body for the 2-body approximation. Not owning.
  not_null<Celestial const*> parent_;
  not_null<Ephemeris<Barycentric>*> const ephemeris_;

  // The past and present trajectory of the body. It ends at |HistoryTime()|
  // unless |*this| was created after |HistoryTime()|, in which case it ends
  // at |current_time_|.  It is advanced with a constant time step.
  std::unique_ptr<DiscreteTrajectory<Barycentric>> history_;

  // The new implementation of history, also encompasses the prolongation.
  DiscreteTrajectory<Barycentric> psychohistory_;
  bool psychohistory_is_authoritative_;

  // A child trajectory of |*history_|. It is forked at |history_->last_time()|
  // and continues until |current_time_|. It is computed with a non-constant
  // timestep, which breaks symplecticity.
  DiscreteTrajectory<Barycentric>* prolongation_ = nullptr;

  // Child trajectory of |*history_|.
  DiscreteTrajectory<Barycentric>* prediction_ = nullptr;

  std::unique_ptr<FlightPlan> flight_plan_;
  bool is_dirty_ = false;

  // We will use union-find algorithms on |Vessel|s.
  not_null<std::unique_ptr<Subset<Vessel>::Node>> const subset_node_;
  friend class Subset<Vessel>::Node;
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

namespace base {

template<>
inline not_null<Subset<ksp_plugin::Vessel>::Node*>
Subset<ksp_plugin::Vessel>::Node::Get(ksp_plugin::Vessel& element) {
  return element.subset_node_.get();
}

}  // namespace base
}  // namespace principia
