#pragma once

#include <vector>
#include <stack>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "ksp_plugin/burn.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using base::not_null;
using geometry::Instant;
using integrators::AdaptiveStepSizeIntegrator;
using physics::DiscreteTrajectory;
using physics::UniqueDiscreteTrajectory;
using physics::Ephemeris;
using quantities::Length;
using quantities::Mass;
using quantities::Speed;

namespace ksp_plugin {

// A stack of |Burn|s that manages a chain of trajectories obtained by executing
// the corresponding |NavigationManœuvre|s.
class FlightPlan {
 public:
  // Creates a |FlightPlan| with no burns whose chain of trajectories is forked
  // from |root| at the latest time prior to |initial_time|, and which starts
  // with the given |initial_mass|.  The trajectories are computed using the
  // given |integrator| in the given |ephemeris|.
  // The given |final_time| and tolerances are used, and may be modified by
  // the mutators |SetFinalTime| and |SetTolerances|.
  FlightPlan(
      not_null<DiscreteTrajectory<Barycentric>*> root,
      Instant const& initial_time,
      Instant const& final_time,
      Mass const& initial_mass,
      not_null<Ephemeris<Barycentric>*> ephemeris,
      AdaptiveStepSizeIntegrator<
          Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance);
  ~FlightPlan() = default;

  int size() const;
  // |index| must be in [0, size()[.
  NavigationManœuvre const& Get(int index);

  // |size()| must be greater than 0.
  void RemoveLast();

  // The following two functions return false and have no effect if the given
  // |burn| would start before |initial_time_| or before the end of the previous
  // burn, or end after |final_time_|.
  bool Append(Burn burn);
  // |size()| must be greater than 0.
  bool ReplaceLast(Burn burn);

  // Returns false and has no effect if |final_time| is before the end of the
  // last manœuvre or before |initial_time_|.
  bool SetFinalTime(Instant const& final_time);

  // Sets the tolerances used to compute the trajectories.  The trajectories are
  // recomputed.
  void SetTolerances(
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance);

 private:
  // Appends |manœuvre| to |manœuvres_|, recomputes the last coast segment until
  // |manœuvre.initial_time()|, and adds a burn and a coast segment.
  // |manœuvre| must fit between |start_of_last_coast()| and |final_time_|.
  void Append(NavigationManœuvre manœuvre);

  // Recomputes all trajectories in |segments_|.
  void RecomputeSegments();

  // Flows the last segment for the duration of |manœuvre| using its intrinsic
  // acceleration.
  void BurnLastSegment(NavigationManœuvre const& manœuvre);
  // Flows the last segment until |final_time| with no intrinsic acceleration.
  void CoastLastSegment(Instant const& final_time);

  // Adds a trajectory to |segment_|, forked at the end of the last one.
  void AddSegment();
  // Forgets the last segment after its fork.
  void ResetLastSegment();

  Instant start_of_last_coast() const;
  Instant start_of_penultimate_coast() const;

  Mass const initial_mass_;
  Instant const initial_time_;
  Instant final_time_;
  Length length_integration_tolerance_;
  Speed speed_integration_tolerance_;
  // Never empty; Starts and ends with a coasting segment; coasting and burning
  // alternate.  Each segment is a fork of the previous one, so the stack
  // structure prevents dangling.
  std::stack<UniqueDiscreteTrajectory<Barycentric>> segments_;
  std::vector<NavigationManœuvre> manœuvres_;
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  AdaptiveStepSizeIntegrator<
      Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator_;
};

}  // namespace ksp_plugin
}  // namespace principia
