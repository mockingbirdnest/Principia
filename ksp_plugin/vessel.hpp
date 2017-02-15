
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

  virtual not_null<Celestial const*> parent() const;
  virtual void set_parent(not_null<Celestial const*> parent);

  virtual void clear_parts();
  virtual void add_part(not_null<Part const*> part);

  virtual DiscreteTrajectory<Barycentric> const& prediction() const;

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

  std::vector<not_null<Part const*>> parts_;

  // The new implementation of history, also encompasses the prolongation.
  DiscreteTrajectory<Barycentric> psychohistory_;
  bool psychohistory_is_authoritative_;

  // Child trajectory of |*history_|.
  DiscreteTrajectory<Barycentric>* prediction_ = nullptr;

  std::unique_ptr<FlightPlan> flight_plan_;

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
