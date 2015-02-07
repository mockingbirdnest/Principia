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
    CHECK_EQ(prolongation_, nullptr);
  }
  return initialized;
}

inline MassiveBody const& Celestial::body() const {
  return *body_;
}

inline bool Celestial::has_parent() const {
  return parent_ != nullptr;
}

inline Celestial const& Celestial::parent() const {
  return *CHECK_NOTNULL(parent_);
}

inline void Celestial::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

inline Trajectory<Barycentric> const& Celestial::history() const {
  CHECK(is_initialized());
  return *history_;
}

inline Trajectory<Barycentric>* Celestial::mutable_history() {
  CHECK(is_initialized());
  return history_.get();
}

inline Trajectory<Barycentric> const& Celestial::prolongation() const {
  CHECK(is_initialized());
  return *prolongation_;
}

inline Trajectory<Barycentric>* Celestial::mutable_prolongation() {
  CHECK(is_initialized());
  return prolongation_;
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
