#pragma once

#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

inline Vessel::Vessel(not_null<Celestial const*> const parent)
    : body_(),
      parent_(parent) {}

inline not_null<MasslessBody const*> Vessel::body() const {
  return &body_;
}

inline bool Vessel::is_synchronized() const {
  bool const synchronized = history_ != nullptr;
  if (synchronized) {
    CHECK(owned_prolongation_ == nullptr);
  }
  return synchronized;
}

inline bool Vessel::is_initialized() const {
  bool const initialized = prolongation_ != nullptr;
  if (!initialized) {
    CHECK(owned_prolongation_ == nullptr);
  }
  return initialized;
}

inline not_null<Celestial const*> Vessel::parent() const {
  return parent_;
}

inline void Vessel::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

inline DiscreteTrajectory<Barycentric> const& Vessel::history() const {
  CHECK(is_synchronized());
  return *history_;
}

inline not_null<DiscreteTrajectory<Barycentric>*> Vessel::mutable_history() {
  CHECK(is_synchronized());
  return history_.get();
}

inline DiscreteTrajectory<Barycentric> const& Vessel::prolongation() const {
  CHECK(is_initialized());
  return *prolongation_;
}

inline not_null<DiscreteTrajectory<Barycentric>*>
Vessel::mutable_prolongation() {
  CHECK(is_initialized());
  return prolongation_;
}

inline std::vector<not_null<DiscreteTrajectory<Barycentric>*>> const&
Vessel::flight_plan() const {
  CHECK(is_initialized());
  return flight_plan_;
}

inline bool Vessel::has_flight_plan() const {
  return !flight_plan_.empty();
}

inline DiscreteTrajectory<Barycentric> const& Vessel::prediction() const {
  CHECK(has_prediction());
  return *prediction_;
}

inline bool Vessel::has_prediction() const {
  return prediction_ != nullptr;
}

inline Vessel::Manœuvres const& Vessel::manœuvres() const {
  return manœuvres_;
}

inline not_null<Vessel::Manœuvres*> Vessel::mutable_manœuvres() {
  return &manœuvres_;
}

inline void Vessel::CreateProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_synchronized());
  CHECK(!is_initialized());
  CHECK(owned_prolongation_ == nullptr);
  owned_prolongation_ = std::make_unique<DiscreteTrajectory<Barycentric>>();
  owned_prolongation_->Append(time, degrees_of_freedom);
  prolongation_ = owned_prolongation_.get();
}

inline void Vessel::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_synchronized());
  history_ = std::make_unique<DiscreteTrajectory<Barycentric>>();
  history_->Append(time, degrees_of_freedom);
  prolongation_ = history_->NewFork(time);
  owned_prolongation_.reset();
}

inline void Vessel::ResetProlongation(Instant const& time) {
  CHECK(is_initialized());
  CHECK(is_synchronized());
  CHECK(owned_prolongation_ == nullptr);
  history_->DeleteFork(&prolongation_);
  prolongation_ = history_->NewFork(time);
}

inline void Vessel::UpdateFlightPlan(
    not_null<Ephemeris<Barycentric>*> ephemeris,
    AdaptiveStepSizeIntegrator<
        Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
    Instant const& last_time,
    Length const& prediction_length_tolerance,
    Speed const& prediction_speed_tolerance,
    Length const& prolongation_length_tolerance,
    Speed const& prolongation_speed_tolerance) {
  if (!is_synchronized()) {
    return;
  }
  DeleteFlightPlan();
  flight_plan_.emplace_back(
      mutable_history()->NewFork(history().last().time()));
  // If prolongation has no additional points this will do nothing (although it
  // might warn).
  flight_plan_.back()->Append(prolongation().last().time(),
                              prolongation().last().degrees_of_freedom());
  for (auto const& manœuvre : manœuvres_) {
    not_null<DiscreteTrajectory<Barycentric>*> coast_trajectory =
        flight_plan_.back();
    ephemeris->FlowWithAdaptiveStep(coast_trajectory,
                                    prediction_length_tolerance,
                                    prediction_speed_tolerance,
                                    integrator,
                                    manœuvre->initial_time());
    flight_plan_.emplace_back(
        coast_trajectory->NewFork(coast_trajectory->last().time()));
    not_null<DiscreteTrajectory<Barycentric>*> burn_trajectory =
        flight_plan_.back();
    burn_trajectory->set_intrinsic_acceleration(manœuvre->acceleration());
    ephemeris->FlowWithAdaptiveStep(burn_trajectory,
                                    prolongation_length_tolerance,
                                    prolongation_speed_tolerance,
                                    integrator,
                                    manœuvre->final_time());
    flight_plan_.emplace_back(
        burn_trajectory->NewFork(burn_trajectory->last().time()));
  }
  ephemeris->FlowWithAdaptiveStep(flight_plan_.back(),
                                  prediction_length_tolerance,
                                  prediction_speed_tolerance,
                                  integrator,
                                  last_time);
}

inline void Vessel::DeleteFlightPlan() {
  if (has_flight_plan()) {
    DiscreteTrajectory<Barycentric>* flight_plan_root = flight_plan_.front();
    flight_plan_.clear();
    history_->DeleteFork(&flight_plan_root);
  }
}

inline void Vessel::UpdatePrediction(
    not_null<Ephemeris<Barycentric>*> ephemeris,
    AdaptiveStepSizeIntegrator<
        Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
    Instant const& last_time,
    Length const& prediction_length_tolerance,
    Speed const& prediction_speed_tolerance) {
  if (!is_synchronized()) {
    return;
  }
  DeletePrediction();
  prediction_ = mutable_history()->NewFork(history().last().time());
  // If prolongation has no additional points this will do nothing (although it
  // might warn).
  prediction_->Append(prolongation().last().time(),
                      prolongation().last().degrees_of_freedom());
  ephemeris->FlowWithAdaptiveStep(prediction_,
                                  prediction_length_tolerance,
                                  prediction_speed_tolerance,
                                  integrator,
                                  last_time);
}

inline void Vessel::DeletePrediction() {
  if (has_prediction()) {
    mutable_history()->DeleteFork(&prediction_);
  }
}

inline void Vessel::WriteToMessage(
    not_null<serialization::Vessel*> const message) const {
  CHECK(is_initialized());
  body_.WriteToMessage(message->mutable_body());
  if (is_synchronized()) {
    history_->WriteToMessage(
        message->mutable_history_and_prolongation()->mutable_history());
    prolongation_->WritePointerToMessage(
        message->mutable_history_and_prolongation()->mutable_prolongation());
  } else {
    owned_prolongation_->WriteToMessage(message->mutable_owned_prolongation());
  }
}

inline std::unique_ptr<Vessel> Vessel::ReadFromMessage(
    serialization::Vessel const& message,
    not_null<Celestial const*> const parent) {
  auto vessel = std::make_unique<Vessel>(parent);
  // NOTE(egg): for now we do not read the |MasslessBody| as it can contain no
  // information.
  if (message.has_history_and_prolongation()) {
    vessel->history_ =
        DiscreteTrajectory<Barycentric>::ReadFromMessage(
            message.history_and_prolongation().history());
    vessel->prolongation_ =
        DiscreteTrajectory<Barycentric>::ReadPointerFromMessage(
            message.history_and_prolongation().prolongation(),
            vessel->history_.get());
  } else if (message.has_owned_prolongation()) {
    vessel->owned_prolongation_ =
        DiscreteTrajectory<Barycentric>::ReadFromMessage(
            message.owned_prolongation());
    vessel->prolongation_ = vessel->owned_prolongation_.get();
  } else {
    LOG(FATAL) << "message does not represent an initialized Vessel";
    base::noreturn();
  }
  return vessel;
}

}  // namespace ksp_plugin
}  // namespace principia
