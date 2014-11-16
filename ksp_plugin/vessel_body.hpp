#pragma once

#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

inline Vessel::Vessel(Celestial const* parent)
    : body_(new Body<Barycentric>(GravitationalParameter())),
      parent_(CHECK_NOTNULL(parent)) {}

inline bool Vessel::has_history() const {
  return history_ != nullptr;
}

inline bool Vessel::has_prolongation() const {
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

inline Trajectory<Barycentric> const& Vessel::prolongation_or_history() const {
    return prolongation_ == nullptr ? *history_ : *prolongation_;
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

inline void Vessel::CreateHistory(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  history_ = std::make_unique<Trajectory<Barycentric>>(*body_);
  history_->Append(time, degrees_of_freedom);
}

inline void Vessel::ResetProlongation(Instant const& time) {
  if (prolongation_ != nullptr) {
    history_->DeleteFork(&prolongation_);
  }
  prolongation_ = history_->Fork(time);
}

}  // namespace ksp_plugin
}  // namespace principia
