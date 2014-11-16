#pragma once

#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

Vessel::Vessel(Celestial const* parent)
    : body_(new Body<Barycentric>(GravitationalParameter())),
      parent_(CHECK_NOTNULL(parent)) {}

bool Vessel::has_history() const {
  return history_ != nullptr;
}

bool Vessel::has_prolongation() const {
  return prolongation_ != nullptr;
}

Celestial const& Vessel::parent() const {
  return *parent_;
}

Trajectory<Barycentric> const& Vessel::history() const {
  return *history_;
}

Trajectory<Barycentric> const& Vessel::prolongation() const {
  return *prolongation_;
}

Trajectory<Barycentric> const& Vessel::prolongation_or_history() const {
    return prolongation_ == nullptr ? *history_ : *prolongation_;
}

Trajectory<Barycentric>* Vessel::mutable_history() {
  return history_.get();
}

Trajectory<Barycentric>* Vessel::mutable_prolongation() {
  return prolongation_;
}

void Vessel::set_parent(Celestial const* parent) {
  parent_ = CHECK_NOTNULL(parent);
}

void Vessel::CreateHistory(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  history_ = std::make_unique<Trajectory<Barycentric>>(*body_);
  history_->Append(time, degrees_of_freedom);
}

void Vessel::ResetProlongation(Instant const& time) {
  if (prolongation_ != nullptr) {
    history_->DeleteFork(&prolongation_);
  }
  prolongation_ = history_->Fork(time);
}

}  // namespace ksp_plugin
}  // namespace principia
