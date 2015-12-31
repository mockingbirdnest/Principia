#pragma once

#include <vector>
#include <stack>

#include "ksp_plugin/fork_handle.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {

using integrators::AdaptiveStepSizeIntegrator;

namespace ksp_plugin {

struct Burn {
  Force thrust;
  SpecificImpulse specific_impulse;
  not_null<std::unique_ptr<NavigationFrame const>> frame;
  Instant initial_time;
  Velocity<Frenet<Navigation>> Δv;
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

  void RemoveLast();

  // The following two functions return false and have no effect if the given
  // |burn| would start before |initial_time_| or before the end of the previous
  // burn, or end after |final_time_|.
  bool Append(Burn burn);
  bool ReplaceLast(Burn burn);

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
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  AdaptiveStepSizeIntegrator<
      Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator_;
};

}  // namespace ksp_plugin
}  // namespace principia
