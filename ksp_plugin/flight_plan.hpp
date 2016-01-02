#pragma once

#include <vector>
#include <stack>

#include "ksp_plugin/burn.hpp"
#include "ksp_plugin/fork_handle.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {

using integrators::AdaptiveStepSizeIntegrator;

namespace ksp_plugin {

using NavigationManœuvre = Manœuvre<Barycentric, Navigation>;

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
  NavigationManœuvre const& Get(int index);

  void RemoveLast();

  // The following two functions return false and have no effect if the given
  // |burn| would start before |initial_time_| or before the end of the previous
  // burn, or end after |final_time_|.
  bool Append(Burn burn);
  bool ReplaceLast(Burn burn);

  // Returns false and has no effect if |final_time| is before the end of the
  // last manœuvre or before |initial_time_|.
  bool SetFinalTime(Instant const& final_time);

  void SetTolerances(
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance);

 private:
  void Append(NavigationManœuvre manœuvre);

  void RecomputeSegments();

  void BurnLastSegment(NavigationManœuvre const& manœuvre);
  void CoastLastSegment(Instant const& final_time);

  void AddSegment();
  void ResetLastSegment();

  // TODO(egg): consider making this a constructor of Manœuvre.
  NavigationManœuvre MakeManœuvre(Burn burn, Mass const& initial_mass);

  Mass const initial_mass_;
  Instant const initial_time_;
  Instant final_time_;
  Length length_integration_tolerance_;
  Speed speed_integration_tolerance_;
  // Starts and ends with a coasting segment; coasting and burning alternate.
  // each segment is a fork of the previous one, so the stack structure prevents
  // dangling.
  std::stack<ForkHandle> segments_;
  std::vector<NavigationManœuvre> manœuvres_;
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  AdaptiveStepSizeIntegrator<
      Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator_;
};

}  // namespace ksp_plugin
}  // namespace principia
