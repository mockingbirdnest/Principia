#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "ksp_plugin/plugin.hpp"

namespace principia {

using base::not_null;

namespace ksp_plugin {

class MockPlugin : public Plugin {
 public:
  MockPlugin();
  MockPlugin(MockPlugin const&) = delete;
  MockPlugin(MockPlugin&&) = delete;
  ~MockPlugin() override = default;

  MOCK_METHOD4(InsertCelestial,
               void(Index const celestial_index,
                    GravitationalParameter const& gravitational_parameter,
                    Index const parent_index,
                    RelativeDegreesOfFreedom<AliceSun> const& from_parent));

  MOCK_METHOD0(EndInitialization,
               void());

  MOCK_CONST_METHOD2(UpdateCelestialHierarchy,
                     void(Index const celestial_index,
                          Index const parent_index));

  MOCK_METHOD2(InsertOrKeepVessel,
               bool(GUID const& vessel_guid, Index const parent_index));

  MOCK_METHOD2(SetVesselStateOffset,
               void(GUID const& vessel_guid,
                    RelativeDegreesOfFreedom<AliceSun> const& from_parent));

  MOCK_METHOD2(AdvanceTime,
               void(Instant const& t, Angle const& planetarium_rotation));

  MOCK_CONST_METHOD1(VesselFromParent,
                     RelativeDegreesOfFreedom<AliceSun>(
                         GUID const& vessel_guid));

  MOCK_CONST_METHOD1(CelestialFromParent,
                     RelativeDegreesOfFreedom<AliceSun>(
                         Index const celestial_index));

  MOCK_CONST_METHOD3(RenderedVesselTrajectory,
                     RenderedTrajectory<World>(
                         GUID const& vessel_guid,
                         not_null<Transforms<
                             Barycentric, Rendering, Barycentric>*> const
                             transforms,
                         Position<World> const& sun_world_position));

  MOCK_METHOD2(RenderedPrediction,
                     RenderedTrajectory<World>(
                         not_null<Transforms<
                             Barycentric, Rendering, Barycentric>*> const
                             transforms,
                         Position<World> const& sun_world_position));

  MOCK_METHOD1(set_predicted_vessel, void(GUID const& vessel_guid));

  MOCK_METHOD0(clear_predicted_vessel, void());

  MOCK_METHOD1(set_prediction_length, void(Time const& t));

  MOCK_METHOD1(set_prediction_step, void(Time const& t));

  // NOTE(phl): gMock 1.7.0 doesn't support returning a std::unique_ptr<>.  So
  // we override the function of the Plugin class with bona fide functions which
  // call mock functions which fill a std::unique_ptr<> instead of returning it.
  not_null<std::unique_ptr<Transforms<Barycentric, Rendering, Barycentric>>>
  NewBodyCentredNonRotatingTransforms(
      Index const reference_body_index) const override;

  not_null<std::unique_ptr<Transforms<Barycentric, Rendering, Barycentric>>>
  NewBarycentricRotatingTransforms(
      Index const primary_index,
      Index const secondary_index) const override;

  MOCK_CONST_METHOD2(FillBodyCentredNonRotatingTransforms,
                     void(Index const reference_body_index,
                          std::unique_ptr<
                              Transforms<Barycentric, Rendering, Barycentric>>*
                                  transforms));

  MOCK_CONST_METHOD3(FillBarycentricRotatingTransforms,
                     void(Index const primary_index,
                          Index const secondary_index,
                          std::unique_ptr<
                              Transforms<Barycentric, Rendering, Barycentric>>*
                                  transforms));

  // NOTE(phl): Another wrapper needed because gMock 1.7.0 wants to copy the
  // vector of unique_ptr<>.
  void AddVesselToNextPhysicsBubble(GUID const& vessel_guid,
                                    std::vector<IdAndOwnedPart> parts) override;

  MOCK_METHOD2(AddVesselToNextPhysicsBubbleConstRef,
               void(GUID const& vessel_guid,
                    std::vector<IdAndOwnedPart> const& parts));

  MOCK_CONST_METHOD0(PhysicsBubbleIsEmpty, bool());

  MOCK_CONST_METHOD1(BubbleDisplacementCorrection,
                     Displacement<World>(
                         Position<World> const& sun_world_position));

  MOCK_CONST_METHOD1(BubbleVelocityCorrection,
                     Velocity<World>(
                         Index const reference_body_index));

  MOCK_CONST_METHOD0(current_time, Instant());

  MOCK_CONST_METHOD1(WriteToMessage,
                     void(not_null<serialization::Plugin*> const message));
};

}  // namespace ksp_plugin
}  // namespace principia
