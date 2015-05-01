#pragma once

#include "ksp_plugin/frames.hpp"
#include "physics/trajectory.hpp"

namespace principia {

using physics::Trajectory;

namespace ksp_plugin {

class MobileInterface {
 public:
  virtual Trajectory<Barycentric> const& history() const = 0;
  virtual not_null<Trajectory<Barycentric>*> mutable_history() = 0;

  virtual Trajectory<Barycentric> const& prolongation() const = 0;
  virtual not_null<Trajectory<Barycentric>*> mutable_prolongation() = 0;

  virtual Trajectory<Barycentric> const& prediction() const = 0;
  virtual Trajectory<Barycentric>* mutable_prediction() = 0;
  virtual bool has_prediction() const = 0;
};

}  // namespace ksp_plugin
}  // namespace principia
