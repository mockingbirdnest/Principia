#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/flight_plan.hpp"

namespace principia {
namespace ksp_plugin {

class MockFlightPlan : public FlightPlan {
 public:
  MockFlightPlan() : FlightPlan() {}

  MOCK_CONST_METHOD0(number_of_manœuvres, int());
  MOCK_CONST_METHOD1(GetManœuvre, NavigationManœuvre const& (int const index));

  MOCK_METHOD0(RemoveLast, void());

  MOCK_CONST_METHOD1(AppendConstRef, bool(Burn const& burn));
  MOCK_CONST_METHOD1(ReplaceLastConstRef, bool(Burn const& burn));

  //TODO(phl):In the .cpp
  bool Append(Burn burn) {
    return AppendConstRef(burn);
  }

  bool ReplaceLast(Burn burn) {
    return ReplaceLastConstRef(burn);
  }

  MOCK_METHOD1(SetFinalTime, bool(Instant const& final_time));

  MOCK_METHOD2(SetTolerances,
               void(Length const& length_integration_tolerance,
                    Speed const& speed_integration_tolerance));

  MOCK_CONST_METHOD0(number_of_segments, int());

  MOCK_CONST_METHOD3(
      GetSegment,
      void(int const index,
           not_null<DiscreteTrajectory<Barycentric>::Iterator*> begin,
           not_null<DiscreteTrajectory<Barycentric>::Iterator*> end));
};

}  // namespace ksp_plugin
}  // namespace principia
