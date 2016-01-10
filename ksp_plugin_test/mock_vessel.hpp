#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

class MockVessel : public Vessel {
 public:
  MockVessel() : Vessel() {}

  MOCK_CONST_METHOD0(body, not_null<MasslessBody const*>());
  MOCK_CONST_METHOD0(is_synchronized, bool());
  MOCK_CONST_METHOD0(is_initialized, bool());

  MOCK_CONST_METHOD0(parent, not_null<Celestial const*>());
  MOCK_METHOD0(set_parent, void(not_null<Celestial const*> const parent));

  MOCK_CONST_METHOD0(history, DiscreteTrajectory<Barycentric> const&());
  MOCK_METHOD0(mutable_history, not_null<DiscreteTrajectory<Barycentric>*>());

  MOCK_CONST_METHOD0(prolongation, DiscreteTrajectory<Barycentric> const&());
  MOCK_METHOD0(mutable_prolongation,
               not_null<DiscreteTrajectory<Barycentric>*>());

  MOCK_CONST_METHOD0(flight_plan, not_null<FlightPlan*>());
  MOCK_CONST_METHOD0(has_flight_plan, bool());

  MOCK_CONST_METHOD0(prediction, DiscreteTrajectory<Barycentric> const&());
  MOCK_CONST_METHOD0(has_prediction, bool());

  MOCK_METHOD2(CreateProlongation,
               void(Instant const& time,
                    DegreesOfFreedom<Barycentric> const& degrees_of_freedom));

  MOCK_METHOD2(CreateHistoryAndForkProlongation,
               void(Instant const& time,
                    DegreesOfFreedom<Barycentric> const& degrees_of_freedom));

  MOCK_METHOD1(ResetProlongation, void(Instant const& time));

  MOCK_METHOD6(CreateFlightPlan,
               void(Instant const& final_time,
                    Mass const& initial_mass,
                    not_null<Ephemeris<Barycentric>*> ephemeris,
                    AdaptiveStepSizeIntegrator<
                        Ephemeris<Barycentric>::NewtonianMotionEquation> const&
                        integrator,
                    Length const& length_integration_tolerance,
                    Speed const& speed_integration_tolerance));

  MOCK_METHOD0(DeleteFlightPlan, void());

  MOCK_METHOD5(UpdatePrediction,
               void(not_null<Ephemeris<Barycentric>*> ephemeris,
                    AdaptiveStepSizeIntegrator<
                        Ephemeris<Barycentric>::NewtonianMotionEquation> const&
                        integrator,
                    Instant const& last_time,
                    Length const& prediction_length_tolerance,
                    Speed const& prediction_speed_tolerance));

  MOCK_METHOD0(DeletePrediction, void());

  MOCK_CONST_METHOD1(WriteToMessage, void(
      not_null<serialization::Vessel*> const message));
};

}  // namespace ksp_plugin
}  // namespace principia
