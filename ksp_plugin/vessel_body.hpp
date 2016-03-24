
#pragma once

#include "ksp_plugin/vessel.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {

using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::si::Kilogram;
using quantities::si::Milli;

namespace ksp_plugin {

inline Vessel::Vessel(not_null<Celestial const*> const parent,
                      not_null<Ephemeris<Barycentric>*> const ephemeris,
                      Ephemeris<Barycentric>::AdaptiveStepParameters const&
                          prolongation_adaptive_step_parameters,
                      Ephemeris<Barycentric>::FixedStepParameters const&
                          history_fixed_step_parameters)
    : body_(),
      parent_(parent),
      ephemeris_(ephemeris),
      prolongation_adaptive_step_parameters_(
          prolongation_adaptive_step_parameters),
      history_fixed_step_parameters_(history_fixed_step_parameters) {}

inline not_null<MasslessBody const*> Vessel::body() const {
  return &body_;
}

inline bool Vessel::is_initialized() const {
  CHECK_EQ(history_ == nullptr, prolongation_ == nullptr);
  return history_ != nullptr;
}

inline not_null<Celestial const*> Vessel::parent() const {
  return parent_;
}

inline void Vessel::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

inline DiscreteTrajectory<Barycentric> const& Vessel::history() const {
  CHECK(is_initialized());
  return *history_;
}

inline DiscreteTrajectory<Barycentric> const& Vessel::prolongation() const {
  CHECK(is_initialized());
  return *prolongation_;
}

inline FlightPlan& Vessel::flight_plan() const {
  CHECK(has_flight_plan());
  return *flight_plan_;
}

inline bool Vessel::has_flight_plan() const {
  return flight_plan_ != nullptr;
}

inline DiscreteTrajectory<Barycentric> const& Vessel::prediction() const {
  CHECK(has_prediction());
  return *prediction_;
}

inline bool Vessel::has_prediction() const {
  return prediction_ != nullptr;
}

inline void Vessel::set_dirty() {
  is_dirty_ = true;
}

inline bool Vessel::is_dirty() const {
  return is_dirty_;
}

inline void Vessel::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_initialized());
  history_ = std::make_unique<DiscreteTrajectory<Barycentric>>();
  history_->Append(time, degrees_of_freedom);
  prolongation_ = history_->NewForkAtLast();
}

inline void Vessel::AdvanceTimeNotInBubble(Instant const& time) {
  CHECK(is_initialized());
  AdvanceHistoryIfNeeded(time);
  FlowProlongation(time);
}

inline void Vessel::AdvanceTimeInBubble(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(is_initialized());
  AdvanceHistoryIfNeeded(time);
  prolongation_->Append(time, degrees_of_freedom);
  is_dirty_ = true;
}

inline void Vessel::ForgetBefore(Instant const& time) {
  CHECK(is_initialized());
  history_->ForgetBefore(time);
}

inline Instant Vessel::ForgettableTime() const {
  CHECK(is_initialized());
  Instant forgettable_time = prolongation_->Fork().time();
  if (flight_plan_ != nullptr) {
    forgettable_time = std::min(forgettable_time,
                                flight_plan_->initial_time());
  }
  if (prediction_ != nullptr) {
    forgettable_time = std::min(forgettable_time,
                                prediction_->Fork().time());
  }
  return forgettable_time;
}

inline void Vessel::CreateFlightPlan(
    Instant const& final_time,
    Mass const& initial_mass,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters) {
  flight_plan_ = std::make_unique<FlightPlan>(
                     history_.get(),
                     /*initial_time=*/history().last().time(),
                     /*final_time=*/final_time,
                     initial_mass,
                     ephemeris_,
                     adaptive_step_parameters);
}

inline void Vessel::DeleteFlightPlan() {
  flight_plan_.reset();
}

inline void Vessel::UpdatePrediction(
    Instant const& last_time,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters) {
  DeletePrediction();
  prediction_ = history_->NewForkAtLast();
  if (history().last().time() != prolongation().last().time()) {
    prediction_->Append(prolongation().last().time(),
                        prolongation().last().degrees_of_freedom());
  }
  ephemeris_->FlowWithAdaptiveStep(
      prediction_,
      Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
      last_time,
      adaptive_step_parameters);
  prediction_last_time_ = last_time;
  prediction_adaptive_step_parameters_ = adaptive_step_parameters;
}

inline void Vessel::DeletePrediction() {
  if (has_prediction()) {
    history_->DeleteFork(&prediction_);
  }
}

inline void Vessel::WriteToMessage(
    not_null<serialization::Vessel*> const message) const {
  CHECK(is_initialized());
  body_.WriteToMessage(message->mutable_body());
  prolongation_adaptive_step_parameters_.WriteToMessage(
      message->mutable_prolongation_adaptive_step_parameters());
  history_fixed_step_parameters_.WriteToMessage(
      message->mutable_history_fixed_step_parameters());
  history_->WriteToMessage(message->mutable_history(),
                           {prolongation_, prediction_});
  if (prediction_last_time_) {
    prediction_last_time_->WriteToMessage(
        message->mutable_prediction_last_time());
  }
  if (prediction_adaptive_step_parameters_) {
    prediction_adaptive_step_parameters_->WriteToMessage(
        message->mutable_prediction_adaptive_step_parameters());
  }
  if (flight_plan_ != nullptr) {
    flight_plan_->WriteToMessage(message->mutable_flight_plan());
  }
  message->set_is_dirty(is_dirty_);
}

inline not_null<std::unique_ptr<Vessel>> Vessel::ReadFromMessage(
    serialization::Vessel const& message,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    not_null<Celestial const*> const parent) {
  // NOTE(egg): for now we do not read the |MasslessBody| as it can contain no
  // information.
  std::unique_ptr<Vessel> vessel;
  bool const is_pre_буняковский = message.has_history_and_prolongation() ||
                                  message.has_owned_prolongation();

  if (is_pre_буняковский) {
    // For the prolongation.
    Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters(
        /*integrator=*/DormandElMikkawyPrince1986RKN434FM<
            Position<Barycentric>>(),
        /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
        /*ength_tolerance=*/1 * Milli(Metre),
        /*speed_tolerance=*/1 * Milli(Metre) / Second);
    // For the history.
    Ephemeris<Barycentric>::FixedStepParameters fixed_step_parameters(
        McLachlanAtela1992Order5Optimal<Position<Barycentric>>(),
        /*step=*/10 * Second);
    vessel = make_not_null_unique<Vessel>(
        parent,
        ephemeris,
        adaptive_step_parameters,
        fixed_step_parameters);

    if (message.has_history_and_prolongation()) {
      vessel->history_ =
          DiscreteTrajectory<Barycentric>::ReadFromMessage(
              message.history_and_prolongation().history(), /*forks=*/{});
      vessel->prolongation_ =
          DiscreteTrajectory<Barycentric>::ReadPointerFromMessage(
              message.history_and_prolongation().prolongation(),
              vessel->history_.get());
      if (message.has_prediction()) {
        vessel->prediction_ =
            DiscreteTrajectory<Barycentric>::ReadPointerFromMessage(
                message.prediction(),
                vessel->history_.get());
      }
      if (message.has_flight_plan()) {
        vessel->flight_plan_ = FlightPlan::ReadFromMessage(
            message.flight_plan(), vessel->history_.get(), ephemeris);
      }
    } else {
      vessel->history_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
                             message.owned_prolongation(), /*forks=*/{});
      vessel->prolongation_ = vessel->history_->NewForkAtLast();
      CHECK(!message.has_prediction());
      CHECK(!message.has_flight_plan());
    }
  } else {
    CHECK(message.has_history() &&
          message.has_history_fixed_step_parameters() &&
          message.has_prolongation_adaptive_step_parameters())
        << message.DebugString();
    vessel = make_not_null_unique<Vessel>(
        parent,
        ephemeris,
        Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
            message.prolongation_adaptive_step_parameters()),
        Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
            message.history_fixed_step_parameters()));
    vessel->history_ =
        DiscreteTrajectory<Barycentric>::ReadFromMessage(
            message.history(),
            {&vessel->prolongation_, &vessel->prediction_});
    if (message.has_prediction_last_time()) {
      vessel->prediction_last_time_ =
          Instant::ReadFromMessage(message.prediction_last_time());
    }
    if (message.has_prediction_adaptive_step_parameters()) {
      vessel->prediction_adaptive_step_parameters_ =
          Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
              message.prediction_adaptive_step_parameters());
    }
    if (vessel->prediction_last_time_ &&
        vessel->prediction_adaptive_step_parameters_) {
      vessel->UpdatePrediction(*vessel->prediction_last_time_,
                               *vessel->prediction_adaptive_step_parameters_);
    }
    if (message.has_flight_plan()) {
      vessel->flight_plan_ = FlightPlan::ReadFromMessage(
          message.flight_plan(), vessel->history_.get(), ephemeris);
    }
    vessel->is_dirty_ = message.is_dirty();
  }
  return std::move(vessel);
}

inline Vessel::Vessel()
    : body_(),
      parent_(testing_utilities::make_not_null<Celestial const*>()),
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()),
      prolongation_adaptive_step_parameters_(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second),
      history_fixed_step_parameters_(
          McLachlanAtela1992Order5Optimal<Position<Barycentric>>(),
          /*step=*/1 * Second) {}

inline void Vessel::AdvanceHistoryIfNeeded(Instant const & time) {
  Instant const& history_last_time = history_->last().time();
  Time const& Δt = history_fixed_step_parameters_.step();

  if (history_last_time + Δt < time) {
    if (is_dirty_) {
      FlowProlongation(history_last_time + Δt);
      history_->Append(history_last_time + Δt,
                       prolongation_->last().degrees_of_freedom());
      is_dirty_ = false;
    }
    FlowHistory(time);
    history_->DeleteFork(&prolongation_);
    prolongation_ = history_->NewForkAtLast();
  }
}

inline void Vessel::FlowHistory(Instant const& time) {
  ephemeris_->FlowWithFixedStep(
      {history_.get()},
      Ephemeris<Barycentric>::kNoIntrinsicAccelerations,
      time,
      history_fixed_step_parameters_);
}

inline void Vessel::FlowProlongation(Instant const& time) {
  Instant const& prolongation_last_time = prolongation_->last().time();
  CHECK_LE(prolongation_last_time, time);
  if (prolongation_last_time == time) {
    return;
  }
  ephemeris_->FlowWithAdaptiveStep(
      prolongation_,
      Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
      time,
      prolongation_adaptive_step_parameters_);
}

}  // namespace ksp_plugin
}  // namespace principia
