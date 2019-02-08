
#pragma once

#include <list>

#include "gmock/gmock.h"
#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

class MockVessel : public Vessel {
 public:
  MockVessel() : Vessel() {}

  MOCK_CONST_METHOD0(body, not_null<MasslessBody const*>());

  MOCK_CONST_METHOD0(parent, not_null<Celestial const*>());
  MOCK_METHOD1(set_parent, void(not_null<Celestial const*> parent));

  MOCK_CONST_METHOD0(psychohistory,
                     std::shared_ptr<DiscreteTrajectory<Barycentric> const>());
  MOCK_CONST_METHOD0(prediction, DiscreteTrajectory<Barycentric> const&());

  MOCK_CONST_METHOD0(flight_plan, FlightPlan&());
  MOCK_CONST_METHOD0(has_flight_plan, bool());

  MOCK_METHOD1(ForgetBefore, void(Instant const& time));

  MOCK_METHOD3(CreateFlightPlan,
               void(Instant const& final_time,
                    Mass const& initial_mass,
                    Ephemeris<Barycentric>::AdaptiveStepParameters const&
                        adaptive_parameters));

  MOCK_METHOD0(DeleteFlightPlan, void());

  MOCK_METHOD1(FlowPrediction, void(Instant const& last_time));

  MOCK_CONST_METHOD0(psychohistory_is_authoritative, bool());

  MOCK_CONST_METHOD1(WriteToMessage,
                     void(not_null<serialization::Vessel*> message));
};

}  // namespace internal_vessel

using internal_vessel::MockVessel;

}  // namespace ksp_plugin
}  // namespace principia
