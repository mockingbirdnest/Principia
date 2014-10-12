#pragma once

#include "ksp_plugin/celestial.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Frame>
Celestial<Frame>::Celestial(
    GravitationalParameter const& gravitational_parameter)
    : body_(new Body<Frame>(gravitational_parameter)) {}

template<typename Frame>
Body<Frame> const& Celestial<Frame>::body() const {
  return *body_;
}

template<typename Frame>
Celestial<Frame> const* Celestial<Frame>::parent() const {
  return parent_;
}

template<typename Frame>
Trajectory<Frame>* Celestial<Frame>::history() const {
  return history_.get();
}

template<typename Frame>
Trajectory<Frame>* Celestial<Frame>::prolongation() const {
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
