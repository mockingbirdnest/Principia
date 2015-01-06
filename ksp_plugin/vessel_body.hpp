#pragma once

#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

inline Vessel::Vessel(Celestial const* parent)
    : body_(),
      parent_(CHECK_NOTNULL(parent)) {}

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

inline void Vessel::set_parent(Celestial const* parent) {
  parent_ = CHECK_NOTNULL(parent);
}

inline void Vessel::CreateProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_synchronized());
  CHECK(!is_initialized());
  CHECK(owned_prolongation_ == nullptr);
  owned_prolongation_ =
      std::make_unique<Trajectory<Barycentric>>(check_not_null(&body_));
  owned_prolongation_->Append(time, degrees_of_freedom);
  prolongation_ = owned_prolongation_.get();
}

inline void Vessel::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_synchronized());
  history_ =
      std::make_unique<Trajectory<Barycentric>>(check_not_null(&body_));
  history_->Append(time, degrees_of_freedom);
  prolongation_ = history_->Fork(time);
  owned_prolongation_.reset();
}

inline void Vessel::ResetProlongation(Instant const& time) {
  CHECK(is_initialized());
  CHECK(is_synchronized());
  CHECK(owned_prolongation_ == nullptr);
  history_->DeleteFork(&prolongation_);
  prolongation_ = history_->Fork(time);
}

}  // namespace ksp_plugin
}  // namespace principia
