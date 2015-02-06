#pragma once

#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

inline Vessel::Vessel(Vessel&& other)  // NOLINT(build/c++11)
    : body_(),
      parent_(std::move(other.parent_)),
      history_(std::move(other.history_)),
      prolongation_(std::move(other.prolongation_)),
      owned_prolongation_(std::move(other.owned_prolongation_)) {}

inline Vessel::Vessel(not_null<Celestial const*> const parent)
    : body_(),
      parent_(parent) {}

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

inline Celestial const& Vessel::parent() const {
  return *parent_;
}

inline Trajectory<Barycentric> const& Vessel::history() const {
  return *history_;
}

inline Trajectory<Barycentric> const& Vessel::prolongation() const {
  return *prolongation_;
}

inline Trajectory<Barycentric>* Vessel::mutable_history() {
  return history_.get();
}

inline Trajectory<Barycentric>* Vessel::mutable_prolongation() {
  return prolongation_;
}

inline void Vessel::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
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

inline Vessel Vessel::ReadFromMessage(serialization::Vessel const& message,
                                      not_null<Celestial const*> const parent) {
  Vessel vessel(parent);
  // NOTE(egg): for now we do not read the |MasslessBody| as it can contain no
  // information.
  if (message.has_history_and_prolongation()) {
    vessel.history_ =
        std::make_unique<Trajectory<Barycentric>>(
            Trajectory<Barycentric>::ReadFromMessage(
                message.history_and_prolongation().history(), &vessel.body_));
    vessel.prolongation_ =
        Trajectory<Barycentric>::ReadPointerFromMessage(
            message.history_and_prolongation().prolongation(),
            vessel.history_.get());
  } else if (message.has_owned_prolongation()) {
    vessel.owned_prolongation_ =
        std::make_unique<Trajectory<Barycentric>>(
            Trajectory<Barycentric>::ReadFromMessage(
                message.owned_prolongation(), &vessel.body_));
    vessel.prolongation_ = vessel.owned_prolongation_.get();
  } else {
    LOG(FATAL) << "message does not represent an initialized Vessel";
    base::noreturn();
  }
  return vessel;
}

}  // namespace ksp_plugin
}  // namespace principia
