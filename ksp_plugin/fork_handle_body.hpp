#pragma once

#include "base/not_null.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace ksp_plugin {

inline ForkHandle::ForkHandle(
    not_null<DiscreteTrajectory<Barycentric>*> trajectory)
    : trajectory_(trajectory) {
  CHECK(!trajectory_->is_root());
}

inline ForkHandle::~ForkHandle() {
  if (trajectory_ != nullptr) {
    trajectory_->parent()->DeleteFork(&trajectory_);
  }
}

inline DiscreteTrajectory<Barycentric>* ForkHandle::operator->() {
  return trajectory_;
}

inline DiscreteTrajectory<Barycentric> const* ForkHandle::operator->() const {
  return trajectory_;
}

inline DiscreteTrajectory<Barycentric>& ForkHandle::operator*() {
  return *trajectory_;
}

inline DiscreteTrajectory<Barycentric> const& ForkHandle::operator*() const {
  return *trajectory_;
}

inline not_null<DiscreteTrajectory<Barycentric>*> ForkHandle::get() {
  return trajectory_;
}

inline not_null<DiscreteTrajectory<Barycentric> const*>
ForkHandle::get() const {
  return trajectory_;
}

}  // namespace ksp_plugin
}  // namespace principia
