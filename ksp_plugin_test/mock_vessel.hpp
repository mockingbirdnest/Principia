
#pragma once

#include <list>

#include "gmock/gmock.h"
#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

class MockVessel : public Vessel {
 public:
  MockVessel() = default;

  MOCK_CONST_METHOD0(body, not_null<MasslessBody const*>());

  MOCK_CONST_METHOD0(parent, not_null<Celestial const*>());
  MOCK_METHOD1(set_parent, void(not_null<Celestial const*> parent));

  MOCK_CONST_METHOD0(psychohistory, DiscreteTrajectory<Barycentric> const&());
  MOCK_CONST_METHOD0(prediction, DiscreteTrajectory<Barycentric> const&());

  MOCK_CONST_METHOD0(flight_plan, FlightPlan&());
  MOCK_CONST_METHOD0(has_flight_plan, bool());

  MOCK_METHOD1(ForgetBefore, void(Instant const& time));

  MOCK_METHOD4(CreateFlightPlan,
               void(Instant const& final_time,
                    Mass const& initial_mass,
                    Ephemeris<Barycentric>::AdaptiveStepParameters const&
                        flight_plan_adaptive_step_parameters,
                    Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
                        const&
                        flight_plan_generalized_adaptive_step_parameters));

  MOCK_METHOD0(DeleteFlightPlan, void());

  MOCK_CONST_METHOD2(WriteToMessage,
                     void(not_null<serialization::Vessel*> message,
                          PileUp::SerializationIndexForPileUp const&
                              serialization_index_for_pile_up));
};

}  // namespace internal_vessel

using internal_vessel::MockVessel;

}  // namespace ksp_plugin
}  // namespace principia
