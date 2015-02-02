#pragma once

#include "ksp_plugin/celestial.hpp"

namespace principia {
namespace ksp_plugin {

inline Celestial::Celestial(not_null<std::unique_ptr<MassiveBody const>> body)
    : body_(std::move(body)) {}

inline MassiveBody const& Celestial::body() const {
  return *body_;
}

inline bool Celestial::has_parent() const {
  return parent_ != nullptr;
}

inline Celestial const& Celestial::parent() const {
  return *CHECK_NOTNULL(parent_);
}

inline Trajectory<Barycentric> const& Celestial::history() const {
  return *history_;
}

inline Trajectory<Barycentric> const& Celestial::prolongation() const {
  return *prolongation_;
}

inline Trajectory<Barycentric>* Celestial::mutable_history() {
  return history_.get();
}

inline Trajectory<Barycentric>* Celestial::mutable_prolongation() {
  return prolongation_;
}

inline void Celestial::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
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

}  // namespace ksp_plugin
}  // namespace principia
