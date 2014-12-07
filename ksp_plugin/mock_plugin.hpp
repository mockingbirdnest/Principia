#pragma once

#include "ksp_plugin/plugin.hpp"

#include "gmock/gmock.h"

namespace principia {
namespace ksp_plugin {

class MockPlugin : public Plugin {
 public:
  MockPlugin() : Plugin(Instant(),
                        Index(0),
                        1 * SIUnit<GravitationalParameter>(),
                        Angle()) {}
  MockPlugin(MockPlugin const&) = delete;
  MockPlugin(MockPlugin&&) = delete;
  ~MockPlugin() override = default;

  MOCK_METHOD5(InsertCelestial,
               void(Index const celestial_index,
                    GravitationalParameter const& gravitational_parameter,
                    Index const parent_index,
                    Displacement<AliceSun> const& from_parent_position,
                    Velocity<AliceSun> const& from_parent_velocity));

  MOCK_METHOD0(EndInitialization,
               void());

  MOCK_CONST_METHOD2(UpdateCelestialHierarchy,
                     void(Index const celestial_index,
                          Index const parent_index));

  MOCK_METHOD2(InsertOrKeepVessel,
               bool(GUID const& vessel_guid, Index const parent_index));

  MOCK_METHOD3(SetVesselStateOffset,
               void(GUID const& vessel_guid,
                    Displacement<AliceSun> const& from_parent_position,
                    Velocity<AliceSun> const& from_parent_velocity));

  MOCK_METHOD2(AdvanceTime,
               void(Instant const& t, Angle const& planetarium_rotation));

  MOCK_CONST_METHOD1(VesselDisplacementFromParent,
                     Displacement<AliceSun>(GUID const& vessel_guid));

  MOCK_CONST_METHOD1(VesselParentRelativeVelocity,
                     Velocity<AliceSun>(GUID const& vessel_guid));

  MOCK_CONST_METHOD1(CelestialDisplacementFromParent,
                     Displacement<AliceSun>(Index const celestial_index));

  MOCK_CONST_METHOD1(CelestialParentRelativeVelocity,
                     Velocity<AliceSun>(Index const celestial_index));

  MOCK_CONST_METHOD3(RenderedVesselTrajectory,
                     RenderedTrajectory<World>(
                         GUID const& vessel_guid,
                         RenderingFrame const& frame,
                         Position<World> const& sun_world_position));

  MOCK_CONST_METHOD1(NewBodyCentredNonRotatingFrame,
                     std::unique_ptr<BodyCentredNonRotatingFrame>(
                         Index const reference_body_index));

  MOCK_CONST_METHOD2(NewBarycentricRotatingFrame,
                     std::unique_ptr<BarycentricRotatingFrame>(
                         Index const primary_index,
                         Index const secondary_index));

  MOCK_CONST_METHOD2(VesselWorldPosition,
                     Position<World>(
                         GUID const& vessel_guid,
                         Position<World> const& parent_world_position));

  MOCK_CONST_METHOD3(VesselWorldVelocity,
                     Velocity<World>(
                         GUID const& vessel_guid,
                         Velocity<World> const& parent_world_velocity,
                         Time const& parent_rotation_period));
};

}  // namespace ksp_plugin
}  // namespace principia
