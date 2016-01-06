#pragma once

#include "base/not_null.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace ksp_plugin {

// Owns a fork.
// TODO(egg): templatize, move to physics, make movable, make NewFork return
// that.
class ForkHandle {
 public:
  ForkHandle(not_null<DiscreteTrajectory<Barycentric>*> trajectory);

  ForkHandle(ForkHandle const&) = delete;
  ForkHandle(ForkHandle&&) = delete;
  ForkHandle& operator=(ForkHandle const&) = delete;
  ForkHandle& operator=(ForkHandle&&) = delete;

  ~ForkHandle();

  DiscreteTrajectory<Barycentric>* operator->();
  DiscreteTrajectory<Barycentric> const* operator->();
  DiscreteTrajectory<Barycentric>& operator*();
  DiscreteTrajectory<Barycentric> const& operator*();

  not_null<DiscreteTrajectory<Barycentric>*> get();
  not_null<DiscreteTrajectory<Barycentric> const*> get();

 private:
  DiscreteTrajectory<Barycentric>* trajectory_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/fork_handle_body.hpp"
