
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "ksp_plugin/plugin.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

class MockPlugin : public Plugin {
 public:
  MockPlugin();
  MockPlugin(MockPlugin const&) = delete;
  MockPlugin(MockPlugin&&) = delete;

  MOCK_METHOD4(
      InsertCelestialAbsoluteCartesian,
      void(Index celestial_index,
           std::optional<Index> const& parent_index,
           serialization::GravityModel::Body const& gravity_model,
           serialization::InitialState::Cartesian::Body const& initial_state));

  MOCK_METHOD0(EndInitialization,
               void());

  MOCK_CONST_METHOD1(HasEncounteredApocalypse,
                     bool(std::string* details));

  MOCK_CONST_METHOD2(UpdateCelestialHierarchy,
                     void(Index celestial_index,
                          Index parent_index));

  MOCK_CONST_METHOD1(CelestialIndexOfBody, Index(MassiveBody const& body));

  MOCK_METHOD5(InsertOrKeepVessel,
               void(GUID const& vessel_guid,
                    std::string const& vessel_name,
                    Index parent_index,
                    bool loaded,
                    bool& inserted));

  MOCK_METHOD4(InsertUnloadedPart,
               void(PartId part_id,
                    std::string const& name,
                    GUID const& vessel_guid,
                    RelativeDegreesOfFreedom<AliceSun> const& from_parent));

  MOCK_METHOD2(AdvanceTime,
               void(Instant const& t, Angle const& planetarium_rotation));

  MOCK_CONST_METHOD1(ForgetAllHistoriesBefore, void(Instant const& t));

  MOCK_CONST_METHOD2(VesselFromParent,
                     RelativeDegreesOfFreedom<AliceSun>(
                         Index parent_index,
                         GUID const& vessel_guid));

  MOCK_CONST_METHOD1(CelestialFromParent,
                     RelativeDegreesOfFreedom<AliceSun>(Index celestial_index));

  MOCK_CONST_METHOD3(CreateFlightPlan,
                     void(GUID const& vessel_guid,
                          Instant const& final_time,
                          Mass const& initial_mass));

  MOCK_CONST_METHOD2(SetPredictionAdaptiveStepParameters,
                     void(GUID const& vessel_guid,
                          Ephemeris<Barycentric>::AdaptiveStepParameters const&
                              prediction_adaptive_step_parameters));

  MOCK_CONST_METHOD1(HasVessel, bool(GUID const& vessel_guid));
  MOCK_CONST_METHOD1(GetVessel, not_null<Vessel*>(GUID const& vessel_guid));

  not_null<std::unique_ptr<Planetarium>> NewPlanetarium(
      Planetarium::Parameters const& parameters,
      Perspective<Navigation, Camera> const& perspective) const override;
  not_null<std::unique_ptr<NavigationFrame>>
  NewBodyCentredNonRotatingNavigationFrame(
      Index reference_body_index) const override;
  not_null<std::unique_ptr<NavigationFrame>>
  NewBarycentricRotatingNavigationFrame(
      Index primary_index,
      Index secondary_index) const override;

  MOCK_CONST_METHOD3(FillPlanetarium,
                     void(Planetarium::Parameters const& parameters,
                          Perspective<Navigation, Camera> const& perspective,
                          std::unique_ptr<Planetarium>* planetarium));
  MOCK_CONST_METHOD2(
      FillBodyCentredNonRotatingNavigationFrame,
      void(Index reference_body_index,
           std::unique_ptr<NavigationFrame>* navigation_frame));
  MOCK_CONST_METHOD3(
      FillBarycentricRotatingNavigationFrame,
      void(Index primary_index,
           Index secondary_index,
           std::unique_ptr<NavigationFrame>* navigation_frame));

  MOCK_CONST_METHOD1(NavballFrameField,
                     std::unique_ptr<FrameField<World, Navball>>(
                         Position<World> const& sun_world_position));

  MOCK_CONST_METHOD1(VesselTangent,
                     Vector<double, World>(GUID const& vessel_guid));

  MOCK_CONST_METHOD1(VesselNormal,
                     Vector<double, World>(GUID const& vessel_guid));

  MOCK_CONST_METHOD1(VesselBinormal,
                     Vector<double, World>(GUID const& vessel_guid));

  MOCK_CONST_METHOD1(VesselVelocity,
                     Velocity<World>(GUID const& vessel_guid));

  MOCK_CONST_METHOD0(CurrentTime, Instant());

  MOCK_CONST_METHOD0(PlanetariumRotation,
                     Rotation<Barycentric, AliceSun> const&());

  MOCK_METHOD0(renderer, Renderer&());
  MOCK_CONST_METHOD0(renderer, Renderer const&());

  MOCK_CONST_METHOD1(WriteToMessage,
                     void(not_null<serialization::Plugin*> message));
};

}  // namespace internal_plugin

using internal_plugin::MockPlugin;

}  // namespace ksp_plugin
}  // namespace principia
