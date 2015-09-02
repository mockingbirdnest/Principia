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

inline Trajectory<Barycentric> const& Vessel::history() const {
  CHECK(is_synchronized());
  return *history_;
}

inline not_null<Trajectory<Barycentric>*> Vessel::mutable_history() {
  CHECK(is_synchronized());
  return history_.get();
}

inline Trajectory<Barycentric> const& Vessel::prolongation() const {
  CHECK(is_initialized());
  return *prolongation_;
}

inline not_null<Trajectory<Barycentric>*> Vessel::mutable_prolongation() {
  CHECK(is_initialized());
  return prolongation_;
}

inline std::vector<not_null<Trajectory<Barycentric>*>> const&
Vessel::predictions() const {
  CHECK(is_initialized());
  return predictions_;
}

inline bool Vessel::has_predictions() const {
  return !predictions_.empty();
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
  owned_prolongation_ = std::make_unique<Trajectory<Barycentric>>(&body_);
  owned_prolongation_->Append(time, degrees_of_freedom);
  prolongation_ = owned_prolongation_.get();
}

inline void Vessel::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_synchronized());
  history_ = std::make_unique<Trajectory<Barycentric>>(&body_);
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

inline void Vessel::UpdatePredictions(
    not_null<Ephemeris<Barycentric>*> ephemeris,
    AdaptiveStepSizeIntegrator<
        Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
    Instant const& last_time,
    Length const& predictions_length_tolerance,
    Speed const& predictions_speed_tolerance,
    Length const& prolongation_length_tolerance,
    Speed const& prolongation_speed_tolerance) {
  DeletePredictions();
  predictions_.emplace_back(
      mutable_prolongation()->NewFork(prolongation().last().time()));
  for (auto const& manœuvre : manœuvres_) {
    not_null<Trajectory<Barycentric>*> parent_trajectory = predictions_.back();
    ephemeris->FlowWithAdaptiveStep(
        parent_trajectory, predictions_length_tolerance,
        predictions_speed_tolerance, integrator, manœuvre->initial_time());
    predictions_.emplace_back(
        mutable_prolongation()->NewFork(parent_trajectory->last().time()));
    not_null<Trajectory<Barycentric>*> child_trajectory = predictions_.back();
    child_trajectory->set_intrinsic_acceleration(manœuvre->acceleration());
    ephemeris->FlowWithAdaptiveStep(
        child_trajectory, prolongation_length_tolerance,
        prolongation_speed_tolerance, integrator, manœuvre->final_time());
    ephemeris->FlowWithAdaptiveStep(
        parent_trajectory, predictions_length_tolerance,
        predictions_speed_tolerance, integrator, last_time);
  }
  ephemeris->FlowWithAdaptiveStep(
      predictions_.back(), predictions_length_tolerance,
      predictions_speed_tolerance, integrator, last_time);
}

inline void Vessel::DeletePredictions() {
  if (has_predictions()) {
    Trajectory<Barycentric>* prediction_root = predictions_.front();
    predictions_.clear();
    prolongation_->DeleteFork(&prediction_root);
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
        Trajectory<Barycentric>::ReadFromMessage(
            message.history_and_prolongation().history(), &vessel->body_);
    vessel->prolongation_ =
        Trajectory<Barycentric>::ReadPointerFromMessage(
            message.history_and_prolongation().prolongation(),
            vessel->history_.get());
  } else if (message.has_owned_prolongation()) {
    vessel->owned_prolongation_ =
        Trajectory<Barycentric>::ReadFromMessage(message.owned_prolongation(),
                                                 &vessel->body_);
    vessel->prolongation_ = vessel->owned_prolongation_.get();
  } else {
    LOG(FATAL) << "message does not represent an initialized Vessel";
    base::noreturn();
  }
  return vessel;
}

}  // namespace ksp_plugin
}  // namespace principia
