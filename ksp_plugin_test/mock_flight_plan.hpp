
#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/flight_plan.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_flight_plan {

class MockFlightPlan : public FlightPlan {
 public:
  MOCK_CONST_METHOD0(initial_time, Instant());
  MOCK_CONST_METHOD0(desired_final_time, Instant());

  MOCK_CONST_METHOD0(number_of_manœuvres, int());
  MOCK_CONST_METHOD1(GetManœuvre, NavigationManœuvre const& (int index));

  MOCK_METHOD0(RemoveLast, void());

  MOCK_METHOD1(Append, bool(NavigationManœuvre::Burn const& burn));
  MOCK_METHOD1(ReplaceLast, bool(NavigationManœuvre::Burn const& burn));

  MOCK_METHOD1(SetDesiredFinalTime, bool(Instant const& final_time));

  MOCK_CONST_METHOD0(adaptive_step_parameters,
                     Ephemeris<Barycentric>::AdaptiveStepParameters const&());
  MOCK_METHOD2(
      SetAdaptiveStepParameters,
      bool(Ephemeris<Barycentric>::AdaptiveStepParameters const&
               adaptive_step_parameters,
           Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
               generalized_adaptive_step_parameters));

  MOCK_METHOD2(SetTolerances,
               void(Length const& length_integration_tolerance,
                    Speed const& speed_integration_tolerance));

  MOCK_CONST_METHOD0(number_of_segments, int());

  MOCK_CONST_METHOD3(GetSegment,
                     void(int index,
                          DiscreteTrajectory<Barycentric>::Iterator& begin,
                          DiscreteTrajectory<Barycentric>::Iterator& end));
};

}  // namespace internal_flight_plan

using internal_flight_plan::MockFlightPlan;

}  // namespace ksp_plugin
}  // namespace principia
