#pragma once

#include "base/not_null.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {

using base::not_null;
using physics::DiscreteTrajectory;

namespace ksp_plugin {

// Owns a fork.
// TODO(egg): templatize, move to physics, make movable, make NewFork return
// that.
class ForkHandle {
 public:
  explicit ForkHandle(not_null<DiscreteTrajectory<Barycentric>*> trajectory);

  ForkHandle(ForkHandle const&) = delete;
  ForkHandle(ForkHandle&&) = delete;
  ForkHandle& operator=(ForkHandle const&) = delete;
  ForkHandle& operator=(ForkHandle&&) = delete;

  ~ForkHandle();

  DiscreteTrajectory<Barycentric>* operator->();
  DiscreteTrajectory<Barycentric> const* operator->() const;
  DiscreteTrajectory<Barycentric>& operator*();
  DiscreteTrajectory<Barycentric> const& operator*() const;

  not_null<DiscreteTrajectory<Barycentric>*> get();
  not_null<DiscreteTrajectory<Barycentric> const*> get() const;

 private:
  DiscreteTrajectory<Barycentric>* trajectory_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/fork_handle_body.hpp"
