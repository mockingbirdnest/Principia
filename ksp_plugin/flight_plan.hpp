#pragma once

#include <vector>
#include <stack>

#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {

using integrators::AdaptiveStepSizeIntegrator;

namespace ksp_plugin {

struct BurnDefinition {
  Force thrust;
  SpecificImpulse specific_impulse;
  not_null<std::unique_ptr<NavigationFrame const>> frame;
  Instant initial_time;
  Velocity<Frenet<Navigation>> Δv;
};

// Owns a fork.
// TODO(egg): templatize, move to physics, make movable, make NewFork return
// that.
class ForkHandle {
 public:
  ForkHandle(not_null<DiscreteTrajectory<Barycentric>*> trajectory)
      : trajectory_(trajectory) {
    CHECK(!trajectory_->is_root());
  }

  ForkHandle(ForkHandle const&) = delete;
  ForkHandle(ForkHandle&&) = delete;
  ForkHandle& operator=(ForkHandle const&) = delete;
  ForkHandle& operator=(ForkHandle&&) = delete;

  ~ForkHandle() {
    if (trajectory_ != nullptr) {
      trajectory_->parent()->DeleteFork(&trajectory_);
    }
  }

  DiscreteTrajectory<Barycentric>* operator->() {
    return trajectory_;
  }

  DiscreteTrajectory<Barycentric> const* operator->() const {
    return trajectory_;
  }

  DiscreteTrajectory<Barycentric>& operator*() {
    return *trajectory_;
  }

  DiscreteTrajectory<Barycentric> const& operator*() const {
    return *trajectory_;
  }

  not_null<DiscreteTrajectory<Barycentric>*> get() {
    return trajectory_;
  }

  not_null<DiscreteTrajectory<Barycentric> const*> get() const {
    return trajectory_;
  }

 private:
  DiscreteTrajectory<Barycentric>* trajectory_;
};

class FlightPlan {
 public:
  FlightPlan(
      not_null<DiscreteTrajectory<Barycentric>*> root,
      Instant const& initial_time,
      Mass const& initial_mass,
      not_null<Ephemeris<Barycentric>*> ephemeris,
      AdaptiveStepSizeIntegrator<
          Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator);
  ~FlightPlan() = default;

  int size() const;
  Manœuvre<Barycentric, Navigation> const& Get(int index);

  // The following two functions return false and have no effect if the given
  // |burn| would start before |initial_time_| or before the end of the previous
  // burn, or end after |final_time_|.
  bool Append(BurnDefinition burn);
  bool ReplaceLast(BurnDefinition burn);

  void RemoveLast();

  // Returns false and has no effect if |final_time| is before the end of the
  // last manœuvre or before |initial_time_|.
  bool set_final_time(Instant const& final_time);

 private:
  Mass const initial_mass_;
  Instant const initial_time_;
  Instant final_time_;
  Length length_integration_tolerance_;
  Speed speed_integration_tolerance_;
  // Starts and ends with a coasting segment; coasting and burning alternate.
  std::stack<ForkHandle> segments_;
  std::vector<Manœuvre<Barycentric, Navigation>> manœuvres_;
  std::stack<BurnDefinition> burns_;
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  AdaptiveStepSizeIntegrator<
      Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator_;
};

}  // namespace ksp_plugin
}  // namespace principia
