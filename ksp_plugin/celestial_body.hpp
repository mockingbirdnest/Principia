#pragma once

#include "ksp_plugin/celestial.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Frame>
template<typename... Args>
Celestial<Frame>::Celestial(Args&&... args)  // NOLINT(build/c++11)
    : body_(new Body<Frame>(std::forward<Args>(args)...)) {}  // NOLINT(build/c++11) NOLINT(whitespace/line_length)

template<typename Frame>
bool Celestial<Frame>::has_parent() const {
  return parent_ != nullptr;
}

template<typename Frame>
Celestial<Frame> const& Celestial<Frame>::parent() const {
  return *CHECK_NOTNULL(parent_);
}

template<typename Frame>
Trajectory<Frame> const& Celestial<Frame>::history() const {
  return *history_;
}

template<typename Frame>
Trajectory<Frame> const& Celestial<Frame>::prolongation() const {
  return *prolongation_;
}

template<typename Frame>
Trajectory<Frame>* Celestial<Frame>::mutable_history() {
  return history_.get();
}

template<typename Frame>
Trajectory<Frame>* Celestial<Frame>::mutable_prolongation() {
  return prolongation_;
}

template<typename Frame>
void Celestial<Frame>::set_parent(Celestial const* parent) {
  parent_ = CHECK_NOTNULL(parent);
}

template<typename Frame>
void Celestial<Frame>::AppendAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  history_ = std::make_unique<Trajectory<Barycentre>>(*body_);
  history_->Append(time, degrees_of_freedom);
  prolongation_ = history_->Fork(time);
}

template<typename Frame>
void Celestial<Frame>::ResetProlongation(Instant const& time) {
  history_->DeleteFork(&prolongation_);
  prolongation_ = history_->Fork(time);
}


}  // namespace ksp_plugin
}  // namespace principia
