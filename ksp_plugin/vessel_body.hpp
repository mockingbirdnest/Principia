#pragma once

#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

inline Vessel::Vessel(Celestial const* parent)
    : body_(new MasslessBody),
      parent_(CHECK_NOTNULL(parent)) {}

inline bool Vessel::synchronized() const {
  return history_ != nullptr;
}

inline bool Vessel::initialized() const {
  return prolongation_ != nullptr;
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
  CHECK(history_ == nullptr);
  CHECK(prolongation_ == nullptr);
  CHECK(owned_prolongation_ == nullptr);
  owned_prolongation_ = std::make_unique<Trajectory<Barycentric>>(*body_);
  owned_prolongation_->Append(time, degrees_of_freedom);
  prolongation_ = owned_prolongation_.get();
}

inline void Vessel::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(history_ == nullptr);
  history_ = std::make_unique<Trajectory<Barycentric>>(*body_);
  history_->Append(time, degrees_of_freedom);
  prolongation_ = history_->Fork(time);
  owned_prolongation_.reset();
}

inline void Vessel::ResetProlongation(Instant const& time) {
  CHECK(history_ != nullptr);
  CHECK(prolongation_ != nullptr);
  CHECK(owned_prolongation_ == nullptr);
  history_->DeleteFork(&prolongation_);
  prolongation_ = history_->Fork(time);
}

}  // namespace ksp_plugin
}  // namespace principia
