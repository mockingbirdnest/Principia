
#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

class MockVessel : public Vessel {
 public:
  MockVessel() : Vessel() {}

  MOCK_CONST_METHOD0(body, not_null<MasslessBody const*>());
  MOCK_CONST_METHOD0(is_initialized, bool());

  MOCK_CONST_METHOD0(parent, not_null<Celestial const*>());
  MOCK_METHOD1(set_parent, void(not_null<Celestial const*> parent));

  MOCK_CONST_METHOD0(history, DiscreteTrajectory<Barycentric> const&());
  MOCK_CONST_METHOD0(prolongation, DiscreteTrajectory<Barycentric> const&());
  MOCK_CONST_METHOD0(prediction, DiscreteTrajectory<Barycentric> const&());

  MOCK_CONST_METHOD0(flight_plan, FlightPlan&());
  MOCK_CONST_METHOD0(has_flight_plan, bool());

  MOCK_METHOD0(set_dirty, void());
  MOCK_CONST_METHOD0(is_dirty, bool());

  MOCK_METHOD2(CreateHistoryAndForkProlongation,
               void(Instant const& time,
                    DegreesOfFreedom<Barycentric> const& degrees_of_freedom));

  MOCK_METHOD1(ForgetBefore, void(Instant const& time));

  MOCK_CONST_METHOD0(ForgettableTime, Instant());

  MOCK_METHOD3(CreateFlightPlan,
               void(Instant const& final_time,
                    Mass const& initial_mass,
                    Ephemeris<Barycentric>::AdaptiveStepParameters const&
                        adaptive_parameters));

  MOCK_METHOD0(DeleteFlightPlan, void());

  MOCK_METHOD1(UpdatePrediction, void(Instant const& last_time));

  MOCK_METHOD0(clear_mass, void());
  MOCK_METHOD1(increment_mass, void(Mass const& mass));
  MOCK_CONST_METHOD0(mass, Mass const&());

  MOCK_METHOD0(clear_intrinsic_force, void());
  MOCK_METHOD1(increment_intrinsic_force,
               void(Vector<Force, Barycentric> const& intrinsic_force));
  MOCK_CONST_METHOD0(intrinsic_force, Vector<Force, Barycentric> const&());

  MOCK_METHOD1(set_containing_pile_up,
               void(IteratorOn<std::list<PileUp>> pile_up));
  MOCK_CONST_METHOD0(
      containing_pile_up,
      std::experimental::optional<IteratorOn<std::list<PileUp>>>());
  MOCK_CONST_METHOD0(is_piled_up, bool());
  MOCK_METHOD0(clear_pile_up, void());

  MOCK_METHOD3(AppendToPsychohistory,
               void(Instant const& time,
                    DegreesOfFreedom<Barycentric> const& degrees_of_freedom,
                    bool authoritative));
  MOCK_CONST_METHOD0(psychohistory, DiscreteTrajectory<Barycentric> const&());
  MOCK_CONST_METHOD0(psychohistory_is_authoritative, bool());

  MOCK_CONST_METHOD1(WriteToMessage,
                     void(not_null<serialization::Vessel*> message));
};

}  // namespace internal_vessel

using internal_vessel::MockVessel;

}  // namespace ksp_plugin
}  // namespace principia
