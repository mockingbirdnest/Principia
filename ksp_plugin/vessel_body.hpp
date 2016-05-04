
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
using quantities::IsFinite;
using quantities::si::Kilogram;
using quantities::si::Milli;

namespace ksp_plugin {

inline Vessel::Vessel(not_null<Celestial const*> const parent,
                      not_null<Ephemeris<Barycentric>*> const ephemeris,
                      Ephemeris<Barycentric>::FixedStepParameters const&
                          history_fixed_step_parameters,
                      Ephemeris<Barycentric>::AdaptiveStepParameters const&
                          prolongation_adaptive_step_parameters,
                      Ephemeris<Barycentric>::AdaptiveStepParameters const&
                          prediction_adaptive_step_parameters)
    : body_(),
      history_fixed_step_parameters_(history_fixed_step_parameters),
      prolongation_adaptive_step_parameters_(
          prolongation_adaptive_step_parameters),
      prediction_adaptive_step_parameters_(
          prediction_adaptive_step_parameters),
      parent_(parent),
      ephemeris_(ephemeris) {}

inline not_null<MasslessBody const*> Vessel::body() const {
  return &body_;
}

inline bool Vessel::is_initialized() const {
  CHECK_EQ(history_ == nullptr, prolongation_ == nullptr);
  CHECK_EQ(history_ == nullptr, prediction_ == nullptr);
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

inline DiscreteTrajectory<Barycentric> const& Vessel::prediction() const {
  CHECK(is_initialized());
  return *prediction_;
}

inline FlightPlan& Vessel::flight_plan() const {
  CHECK(has_flight_plan());
  return *flight_plan_;
}

inline bool Vessel::has_flight_plan() const {
  return flight_plan_ != nullptr;
}

inline void Vessel::set_dirty() {
  is_dirty_ = true;
}

inline bool Vessel::is_dirty() const {
  return is_dirty_;
}

inline void Vessel::set_prediction_adaptive_step_parameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        prediction_adaptive_step_parameters) {
  prediction_adaptive_step_parameters_ = prediction_adaptive_step_parameters;
}

inline Ephemeris<Barycentric>::AdaptiveStepParameters const&
Vessel::prediction_adaptive_step_parameters() const {
  return prediction_adaptive_step_parameters_;
}

inline void Vessel::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_initialized());
  history_ = std::make_unique<DiscreteTrajectory<Barycentric>>();
  history_->Append(time, degrees_of_freedom);
  prolongation_ = history_->NewForkAtLast();
  prediction_ = history_->NewForkAtLast();
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
  if (prediction_->Fork().time() < time) {
    history_->DeleteFork(&prediction_);
    prediction_ = history_->NewForkAtLast();
  }
  if (flight_plan_ != nullptr) {
    flight_plan_->ForgetBefore(time, [this]() { flight_plan_.reset(); });
  }
  history_->ForgetBefore(time);
}

inline void Vessel::CreateFlightPlan(
    Instant const& final_time,
    Mass const& initial_mass,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        flight_plan_adaptive_step_parameters) {
  auto const history_last = history().last();
  flight_plan_ = std::make_unique<FlightPlan>(
      initial_mass,
      /*initial_time=*/history_last.time(),
      /*initial_degrees_of_freedom=*/history_last.degrees_of_freedom(),
      final_time,
      ephemeris_,
      flight_plan_adaptive_step_parameters);
}

inline void Vessel::DeleteFlightPlan() {
  flight_plan_.reset();
}

inline void Vessel::UpdatePrediction(Instant const& last_time) {
  CHECK(is_initialized());
  history_->DeleteFork(&prediction_);
  prediction_ = history_->NewForkAtLast();
  auto const prolongation_last = prolongation_->last();
  if (history_->last().time() != prolongation_last.time()) {
    prediction_->Append(prolongation_last.time(),
                        prolongation_last.degrees_of_freedom());
  }
  FlowPrediction(last_time);
}

inline void Vessel::WriteToMessage(
    not_null<serialization::Vessel*> const message) const {
  CHECK(is_initialized());
  body_.WriteToMessage(message->mutable_body());
  prolongation_adaptive_step_parameters_.WriteToMessage(
      message->mutable_prolongation_adaptive_step_parameters());
  history_fixed_step_parameters_.WriteToMessage(
      message->mutable_history_fixed_step_parameters());
  history_->WriteToMessage(message->mutable_history(), {prolongation_});
  prediction_->Fork().time().WriteToMessage(
      message->mutable_prediction_fork_time());
  prediction_->last().time().WriteToMessage(
      message->mutable_prediction_last_time());
  prediction_adaptive_step_parameters_.WriteToMessage(
      message->mutable_prediction_adaptive_step_parameters());
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
    vessel = make_not_null_unique<Vessel>(parent,
                                          ephemeris,
                                          DefaultHistoryParameters(),
                                          DefaultProlongationParameters(),
                                          DefaultPredictionParameters());

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
    if (vessel->prediction_ == nullptr) {
      vessel->prediction_ = vessel->history_->NewForkAtLast();
    }
  } else {
    CHECK(message.has_history() &&
          message.has_history_fixed_step_parameters() &&
          message.has_prolongation_adaptive_step_parameters() &&
          message.has_prediction_fork_time() &&
          message.has_prediction_last_time() &&
          message.has_prediction_adaptive_step_parameters())
        << message.DebugString();
    vessel = make_not_null_unique<Vessel>(
        parent,
        ephemeris,
        Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
            message.history_fixed_step_parameters()),
        Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
            message.prolongation_adaptive_step_parameters()),
        Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
            message.prediction_adaptive_step_parameters()));
    vessel->history_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(), {&vessel->prolongation_});
    vessel->prediction_ = vessel->history_->NewForkWithoutCopy(
        Instant::ReadFromMessage(message.prediction_fork_time()));
    vessel->FlowPrediction(
        Instant::ReadFromMessage(message.prediction_last_time()));
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
      history_fixed_step_parameters_(DefaultHistoryParameters()),
      prolongation_adaptive_step_parameters_(DefaultProlongationParameters()),
      prediction_adaptive_step_parameters_(DefaultPredictionParameters()),
      parent_(testing_utilities::make_not_null<Celestial const*>()),
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()) {}

inline void Vessel::AdvanceHistoryIfNeeded(Instant const& time) {
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
      prolongation_adaptive_step_parameters_,
      Ephemeris<Barycentric>::unlimited_max_ephemeris_steps);
}

inline void Vessel::FlowPrediction(Instant const& time) {
  if (time > prediction_->last().time()) {
    bool finite_time = IsFinite(time - prediction_->last().time());
    Instant const t = finite_time ? time : ephemeris_->t_max();
    // This will not prolong the ephemeris if |time| is infinite (but it may do
    // no if it is finite).
    bool reached_t = ephemeris_->FlowWithAdaptiveStep(
        prediction_,
        Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
        t,
        prediction_adaptive_step_parameters_,
        FlightPlan::max_ephemeris_steps_per_frame);
    if (!finite_time && reached_t) {
      // This will prolong the ephemeris by |max_ephemeris_steps_per_frame|.
      ephemeris_->FlowWithAdaptiveStep(
        prediction_,
        Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
        time,
        prediction_adaptive_step_parameters_,
        FlightPlan::max_ephemeris_steps_per_frame);
    }
  }
}

inline Ephemeris<Barycentric>::FixedStepParameters DefaultHistoryParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
             McLachlanAtela1992Order5Optimal<Position<Barycentric>>(),
             /*step=*/10 * Second);
}

inline Ephemeris<Barycentric>::AdaptiveStepParameters
DefaultProlongationParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
             DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
             /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
             /*length_integration_tolerance=*/1 * Milli(Metre),
             /*speed_integration_tolerance=*/1 * Milli(Metre) / Second);
}

inline Ephemeris<Barycentric>::AdaptiveStepParameters
DefaultPredictionParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
             DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
             /*max_steps=*/1000,
             /*length_integration_tolerance=*/1 * Metre,
             /*speed_integration_tolerance=*/1 * Metre / Second);
}

}  // namespace ksp_plugin
}  // namespace principia
