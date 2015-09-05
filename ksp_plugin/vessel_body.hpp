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

inline Trajectory<Barycentric> const& Vessel::prediction() const {
  CHECK(is_initialized());
  return *CHECK_NOTNULL(prediction_);
}

inline Trajectory<Barycentric>* Vessel::mutable_prediction() {
  CHECK(is_initialized());
  return prediction_;
}

inline bool Vessel::has_prediction() const {
  return prediction_ != nullptr;
}

inline void Vessel::CreateProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_synchronized());
  CHECK(!is_initialized());
  CHECK(owned_prolongation_ == nullptr);
  owned_prolongation_ = std::make_unique<Trajectory<Barycentric>>();
  owned_prolongation_->Append(time, degrees_of_freedom);
  prolongation_ = owned_prolongation_.get();
}

inline void Vessel::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_synchronized());
  history_ = std::make_unique<Trajectory<Barycentric>>();
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

inline void Vessel::ForkPrediction() {
  CHECK(prediction_ == nullptr);
  prediction_ = mutable_prolongation()->NewFork(prolongation().last().time());
}

inline void Vessel::DeletePrediction() {
  prolongation_->DeleteFork(&prediction_);
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
            message.history_and_prolongation().history());
    vessel->prolongation_ =
        Trajectory<Barycentric>::ReadPointerFromMessage(
            message.history_and_prolongation().prolongation(),
            vessel->history_.get());
  } else if (message.has_owned_prolongation()) {
    vessel->owned_prolongation_ =
        Trajectory<Barycentric>::ReadFromMessage(message.owned_prolongation());
    vessel->prolongation_ = vessel->owned_prolongation_.get();
  } else {
    LOG(FATAL) << "message does not represent an initialized Vessel";
    base::noreturn();
  }
  return vessel;
}

}  // namespace ksp_plugin
}  // namespace principia
