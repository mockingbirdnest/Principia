
#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/flight_plan.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_flight_plan {

class MockFlightPlan : public FlightPlan {
 public:
  MOCK_METHOD(Instant, initial_time, (), (const, override));
  MOCK_METHOD(Instant, desired_final_time, (), (const, override));

  MOCK_METHOD(int, number_of_manœuvres, (), (const, override));
  MOCK_METHOD(NavigationManœuvre const&,
              GetManœuvre,
              (int index),
              (const, override));

  MOCK_METHOD(absl::Status,
              Insert,
              (NavigationManœuvre::Burn const& burn, int index),
              (override));
  MOCK_METHOD(absl::Status, Remove, (int index), (override));
  MOCK_METHOD(absl::Status,
              Replace,
              (NavigationManœuvre::Burn const& burn, int index),
              (override));

  MOCK_METHOD(absl::Status,
              SetDesiredFinalTime,
              (Instant const& final_time),
              (override));

  MOCK_METHOD(Ephemeris<Barycentric>::AdaptiveStepParameters const&,
              adaptive_step_parameters,
              (),
              (const, override));
  MOCK_METHOD(Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&,
              generalized_adaptive_step_parameters,
              (),
              (const, override));
  MOCK_METHOD(absl::Status,
              SetAdaptiveStepParameters,
              (Ephemeris<Barycentric>::AdaptiveStepParameters const&
                   adaptive_step_parameters,
               Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
                   generalized_adaptive_step_parameters),
              (override));

  MOCK_METHOD(int, number_of_segments, (), (const, override));

  MOCK_METHOD(void,
              GetSegment,
              (int index,
               DiscreteTrajectory<Barycentric>::Iterator& begin,
               DiscreteTrajectory<Barycentric>::Iterator& end),
              (const, override));
};

}  // namespace internal_flight_plan

using internal_flight_plan::MockFlightPlan;

}  // namespace ksp_plugin
}  // namespace principia
