﻿
#pragma once

#include <vector>

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
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using base::not_null;
using geometry::Instant;
using integrators::AdaptiveStepSizeIntegrator;
using physics::DiscreteTrajectory;
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
      not_null<DiscreteTrajectory<Barycentric>*> const root,
      Instant const& initial_time,
      Instant const& final_time,
      Mass const& initial_mass,
      not_null<Ephemeris<Barycentric>*> const ephemeris,
      AdaptiveStepSizeIntegrator<
          Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance);
  virtual ~FlightPlan() = default;

  virtual Instant initial_time() const;
  virtual Instant final_time() const;

  virtual int number_of_manœuvres() const;
  // |index| must be in [0, number_of_manœuvres()[.
  virtual NavigationManœuvre const& GetManœuvre(int const index) const;

  // |size()| must be greater than 0.
  virtual void RemoveLast();

  // The following two functions return false and have no effect if the given
  // |burn| would start before |initial_time_| or before the end of the previous
  // burn, or end after |final_time_|, or if the integration of the coasting
  // phase times out or is singular before the burn.
  virtual bool Append(Burn burn);
  // |size()| must be greater than 0.
  virtual bool ReplaceLast(Burn burn);

  // Returns false and has no effect if |final_time| is before the end of the
  // last manœuvre or before |initial_time_|.
  virtual bool SetFinalTime(Instant const& final_time);

  // Sets the tolerances used to compute the trajectories.  The trajectories are
  // recomputed.
  virtual void SetTolerances(
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance);

  // Returns the number of trajectory segments in this object.
  virtual int number_of_segments() const;

  // |index| must be in [0, number_of_segments()[.  Sets the iterators to denote
  // the given trajectory segment.
  virtual void GetSegment(
      int const index,
      not_null<DiscreteTrajectory<Barycentric>::Iterator*> begin,
      not_null<DiscreteTrajectory<Barycentric>::Iterator*> end) const;

  void WriteToMessage(not_null<serialization::FlightPlan*> const message) const;
  static std::unique_ptr<FlightPlan> ReadFromMessage(
      serialization::FlightPlan const& message,
      not_null<DiscreteTrajectory<Barycentric>*> const root,
      not_null<Ephemeris<Barycentric>*> const ephemeris);

 protected:
  // For mocking.
  FlightPlan();

 private:
  // Appends |manœuvre| to |manœuvres_|, adds a burn and a coast segment.
  // |manœuvre| must fit between |start_of_last_coast()| and |final_time_|,
  // the last coast segment must end at |manœuvre.initial_time()|.
  void Append(NavigationManœuvre manœuvre);

  // Recomputes all trajectories in |segments_|.
  void RecomputeSegments();

  // Flows the last segment for the duration of |manœuvre| using its intrinsic
  // acceleration.
  void BurnLastSegment(NavigationManœuvre const& manœuvre);
  // Flows the last segment until |final_time| with no intrinsic acceleration.
  void CoastLastSegment(Instant const& final_time);

  // Replaces the last segment with |segment|.  |segment| must be forked from
  // the same trajectory as the last segment, and at the same time.  |segment|
  // must not be anomalous.
  void ReplaceLastSegment(
      not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> segment);

  // Adds a trajectory to |segments_|, forked at the end of the last one.
  void AddSegment();
  // Forgets the last segment after its fork.
  void ResetLastSegment();

  // Deletes the last segment and removes it from |segments_|.
  void PopLastSegment();

  // If the integration of a coast from the fork of |coast| until
  // |manœuvre.initial_time()| reaches the end, returns the integrated
  // trajectory.  Otherwise, returns null.
  std::unique_ptr<DiscreteTrajectory<Barycentric>>
  CoastIfReachesManœuvreInitialTime(
      DiscreteTrajectory<Barycentric> const& coast,
      NavigationManœuvre const& manœuvre) const;

  Instant start_of_last_coast() const;
  Instant start_of_penultimate_coast() const;

  DiscreteTrajectory<Barycentric> const& last_coast() const;
  DiscreteTrajectory<Barycentric> const& penultimate_coast() const;

  Mass const initial_mass_;
  Instant initial_time_;
  Instant final_time_;
  Length length_integration_tolerance_;
  Speed speed_integration_tolerance_;
  // Never empty; Starts and ends with a coasting segment; coasting and burning
  // alternate.  This simulates a stack.  Every segment is a nonempty root.
  // For every segment |s| following a segment |p|,
  // |s->begin().time() == p->last().time()|.
  std::vector<not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>>>
      segments_;
  std::vector<NavigationManœuvre> manœuvres_;
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  AdaptiveStepSizeIntegrator<
      Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator_;
  // The last |anomalous_segments_| of |segments_| are anomalous, i.e. they
  // either end prematurely or follow an anomalous segment; in the latter case
  // they are empty.
  // The contract of |Append| and |ReplaceLast| implies that
  // |anomalous_segments_| is at most 2: the penultimate coast is never
  // anomalous.
  int anomalous_segments_ = 0;
};

}  // namespace ksp_plugin
}  // namespace principia
