
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

  MOCK_METHOD2(InsertSun,
               void(Index const celestial_index,
                    GravitationalParameter const& gravitational_parameter));

  void DirectlyInsertCelestial(
      Index const celestial_index,
      Index const* const parent_index,
      DegreesOfFreedom<Barycentric> const& initial_state,
      std::unique_ptr<MassiveBody> body) override;

  MOCK_METHOD4(DirectlyInsertCelestialConstRef,
               void(Index const celestial_index,
                    Index const* const parent_index,
                    DegreesOfFreedom<Barycentric> const& initial_state,
                    std::unique_ptr<MassiveBody> const& body));

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

  MOCK_CONST_METHOD1(ForgetAllHistoriesBefore, void(Instant const& t));

  MOCK_CONST_METHOD1(VesselFromParent,
                     RelativeDegreesOfFreedom<AliceSun>(
                         GUID const& vessel_guid));

  MOCK_CONST_METHOD1(CelestialFromParent,
                     RelativeDegreesOfFreedom<AliceSun>(
                         Index const celestial_index));

  MOCK_CONST_METHOD3(CreateFlightPlan,
                     void(GUID const& vessel_guid,
                          Instant const& final_time,
                          Mass const& initial_mass));

  MOCK_CONST_METHOD2(
      RenderedVesselTrajectory,
      RenderedTrajectory<World>(GUID const& vessel_guid,
                                Position<World> const& sun_world_position));

  MOCK_CONST_METHOD2(
      RenderedPrediction,
      RenderedTrajectory<World>(GUID const& vessel_guid,
                                Position<World> const& sun_world_position));

  MOCK_CONST_METHOD3(
      RenderedTrajectoryFromIterators,
      RenderedTrajectory<World>(
          DiscreteTrajectory<Barycentric>::Iterator const& begin,
          DiscreteTrajectory<Barycentric>::Iterator const& end,
          Position<World> const& sun_world_position));

  MOCK_METHOD1(SetPredictionLength, void(Time const& t));

  MOCK_METHOD1(SetPredictionLengthTolerance, void(Length const& t));

  MOCK_METHOD1(SetPredictionSpeedTolerance, void(Speed const& t));

  MOCK_CONST_METHOD1(HasVessel, bool(GUID const& vessel_guid));
  MOCK_CONST_METHOD1(GetVessel, not_null<Vessel*>(GUID const& vessel_guid));

  // NOTE(phl): gMock 1.7.0 doesn't support returning a std::unique_ptr<>.  So
  // we override the function of the Plugin class with bona fide functions which
  // call mock functions which fill a std::unique_ptr<> instead of returning it.
  not_null<std::unique_ptr<NavigationFrame>>
  NewBodyCentredNonRotatingNavigationFrame(
      Index const reference_body_index) const override;

  not_null<std::unique_ptr<NavigationFrame>>
  NewBarycentricRotatingNavigationFrame(
      Index const primary_index,
      Index const secondary_index) const override;

  MOCK_CONST_METHOD2(
      FillBodyCentredNonRotatingNavigationFrame,
      void(Index const reference_body_index,
           std::unique_ptr<NavigationFrame>* navigation_frame));

  MOCK_CONST_METHOD3(
      FillBarycentricRotatingNavigationFrame,
      void(Index const primary_index,
           Index const secondary_index,
           std::unique_ptr<NavigationFrame>* navigation_frame));

  // NOTE(phl): Needed because gMock 1.7.0 wants to copy the unique_ptr<>.
  void SetPlottingFrame(
    not_null<std::unique_ptr<NavigationFrame>> plotting_frame) override;

  MOCK_METHOD1(SetPlottingFrameConstRef,
               void(NavigationFrame const& plotting_frame));

  MOCK_CONST_METHOD0(GetPlottingFrame,
                     not_null<NavigationFrame const*>());

  MOCK_CONST_METHOD3(
      PlotBarycentricPosition,
      Position<World>(Instant const& t,
                      Position<Barycentric> const& position,
                      Position<World> const& sun_world_position));

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

  MOCK_CONST_METHOD1(Navball,
                     FrameField<World>(
                         Position<World> const& sun_world_position));

  MOCK_CONST_METHOD1(VesselTangent,
                     Vector<double, World>(GUID const& vessel_guid));

  MOCK_CONST_METHOD1(VesselNormal,
                     Vector<double, World>(GUID const& vessel_guid));

  MOCK_CONST_METHOD1(VesselBinormal,
                     Vector<double, World>(GUID const& vessel_guid));

  MOCK_CONST_METHOD0(BarycentricToWorldSun,
                     OrthogonalMap<Barycentric, WorldSun>());

  MOCK_CONST_METHOD0(CurrentTime, Instant());

  MOCK_CONST_METHOD1(WriteToMessage,
                     void(not_null<serialization::Plugin*> const message));
};

}  // namespace ksp_plugin
}  // namespace principia
