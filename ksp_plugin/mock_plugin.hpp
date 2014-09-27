#include "ksp_plugin/plugin.hpp"

#include "gmock/gmock.h"

namespace principia {
namespace ksp_plugin {

class MockPlugin : public Plugin {
 public:
  MockPlugin() = delete;
  MockPlugin(MockPlugin const&) = delete;
  MockPlugin(MockPlugin&&) = delete;
  ~MockPlugin() override = default;

  MockPlugin(Instant const& initial_time,
             Index const sun_index,
             GravitationalParameter const& sun_gravitational_parameter,
             Angle const& planetarium_rotation);

  MOCK_METHOD5(InsertCelestial,
               void(Index const celestial_index,
                    GravitationalParameter const& gravitational_parameter,
                    Index const parent_index,
                    Displacement<AliceSun> const& from_parent_position,
                    Velocity<AliceSun> const& from_parent_velocity));

  MOCK_CONST_METHOD2(UpdateCelestialHierarchy,
                     void(Index const celestial_index,
                          Index const parent_index));

  MOCK_METHOD2(InsertOrKeepVessel,
               bool(GUID const& vessel_guid, Index const parent_index));

  MOCK_CONST_METHOD3(SetVesselStateOffset,
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
                     Displacement<AliceSun> (Index const celestial_index));

  MOCK_CONST_METHOD1(CelestialParentRelativeVelocity,
                     Velocity<AliceSun> (Index const celestial_index));
};

}  // namespace ksp_plugin
}  // namespace principia
