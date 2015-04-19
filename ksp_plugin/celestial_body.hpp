#pragma once

#include "ksp_plugin/celestial.hpp"

namespace principia {
namespace ksp_plugin {

inline Celestial::Celestial(not_null<std::unique_ptr<MassiveBody const>> body)
    : body_(std::move(body)) {}

inline bool Celestial::is_initialized() const {
  bool const initialized = history_ != nullptr;
  if (initialized) {
    CHECK_NOTNULL(prolongation_);
  } else {
    CHECK(prolongation_ == nullptr);
  }
  return initialized;
}

inline MassiveBody const& Celestial::body() const {
  return *body_;
}

inline bool Celestial::has_parent() const {
  return parent_ != nullptr;
}

inline Celestial const* Celestial::parent() const {
  return parent_;
}

inline void Celestial::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

inline Trajectory<Barycentric> const& Celestial::history() const {
  CHECK(is_initialized());
  return *history_;
}

inline not_null<Trajectory<Barycentric>*> Celestial::mutable_history() {
  CHECK(is_initialized());
  return history_.get();
}

inline Trajectory<Barycentric> const& Celestial::prolongation() const {
  CHECK(is_initialized());
  return *prolongation_;
}

inline not_null<Trajectory<Barycentric>*> Celestial::mutable_prolongation() {
  CHECK(is_initialized());
  return prolongation_;
}

inline Trajectory<Barycentric> const& Celestial::prediction() const {
  CHECK(is_initialized());
  return *CHECK_NOTNULL(prediction_);
}

inline Trajectory<Barycentric>* Celestial::mutable_prediction() {
  CHECK(is_initialized());
  return prediction_;
}

inline void Celestial::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  history_ = std::make_unique<Trajectory<Barycentric>>(body_.get());
  history_->Append(time, degrees_of_freedom);
  prolongation_ = history_->NewFork(time);
}

inline void Celestial::ResetProlongation(Instant const& time) {
  history_->DeleteFork(&prolongation_);
  prolongation_ = history_->NewFork(time);
}

inline void Celestial::ForkPrediction() {
  CHECK(prediction_ == nullptr);
  prediction_ = mutable_prolongation()->NewFork(prolongation().last().time());
}

inline void Celestial::DeletePrediction() {
  prolongation_->DeleteFork(&prediction_);
}

inline void Celestial::WriteToMessage(
    not_null<serialization::Celestial*> const message) const {
  CHECK(is_initialized());
  body_->WriteToMessage(message->mutable_body());
  history_->WriteToMessage(
      message->mutable_history_and_prolongation()->mutable_history());
  prolongation_->WritePointerToMessage(
      message->mutable_history_and_prolongation()->mutable_prolongation());
}

inline std::unique_ptr<Celestial> Celestial::ReadFromMessage(
    serialization::Celestial const& message) {
  auto celestial =
      std::make_unique<Celestial>(MassiveBody::ReadFromMessage(message.body()));
  celestial->history_ =
      Trajectory<Barycentric>::ReadFromMessage(
          message.history_and_prolongation().history(),
          celestial->body_.get());
  celestial->prolongation_ =
      Trajectory<Barycentric>::ReadPointerFromMessage(
          message.history_and_prolongation().prolongation(),
          celestial->history_.get());
  return celestial;
}

}  // namespace ksp_plugin
}  // namespace principia
