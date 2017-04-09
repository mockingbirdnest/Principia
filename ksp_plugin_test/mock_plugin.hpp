
#pragma once

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

  void InsertCelestialAbsoluteCartesian(
      Index celestial_index,
      std::experimental::optional<Index> const& parent_index,
      DegreesOfFreedom<Barycentric> const& initial_state,
      base::not_null<std::unique_ptr<RotatingBody<Barycentric> const>> body)
      override;

  MOCK_METHOD4(
      InsertCelestialAbsoluteCartesianConstRef,
      void(Index celestial_index,
           std::experimental::optional<Index> const& parent_index,
           DegreesOfFreedom<Barycentric> const& initial_state,
           base::not_null<
               std::unique_ptr<RotatingBody<Barycentric> const>> const& body));

  MOCK_METHOD0(EndInitialization,
               void());

  MOCK_CONST_METHOD1(HasEncounteredApocalypse,
                     bool(std::string* details));

  MOCK_CONST_METHOD2(UpdateCelestialHierarchy,
                     void(Index celestial_index,
                          Index parent_index));

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

  // NOTE(phl): gMock 1.7.0 doesn't support returning a std::unique_ptr<>.  So
  // we override the function of the Plugin class with bona fide functions which
  // call mock functions which fill a std::unique_ptr<> instead of returning it.
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> RenderedVesselTrajectory(
      GUID const& vessel_guid,
      Position<World> const& sun_world_position) const override;
  MOCK_CONST_METHOD3(FillRenderedVesselTrajectory,
                     void(GUID const& vessel_guid,
                          Position<World> const& sun_world_position,
                          std::unique_ptr<DiscreteTrajectory<World>>*
                              rendered_vessel_trajectory));

  not_null<std::unique_ptr<DiscreteTrajectory<World>>> RenderedPrediction(
      GUID const& vessel_guid,
      Position<World> const& sun_world_position) const override;
  MOCK_CONST_METHOD3(
      FillRenderedPrediction,
      void(GUID const& vessel_guid,
           Position<World> const& sun_world_position,
           std::unique_ptr<DiscreteTrajectory<World>>* rendered_prediction));

  not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderedTrajectoryFromIterators(
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end,
      Position<World> const& sun_world_position) const;
  MOCK_CONST_METHOD4(
      FillRenderedTrajectoryFromIterators,
      void(DiscreteTrajectory<Barycentric>::Iterator const& begin,
           DiscreteTrajectory<Barycentric>::Iterator const& end,
           Position<World> const& sun_world_position,
           std::unique_ptr<DiscreteTrajectory<World>>*
               rendered_trajectory_from_iterators));

  MOCK_METHOD1(SetPredictionLength, void(Time const& t));

  MOCK_METHOD1(SetPredictionAdaptiveStepParameters,
               void(Ephemeris<Barycentric>::AdaptiveStepParameters const&
                        prediction_adaptive_step_parameters));

  MOCK_CONST_METHOD1(HasVessel, bool(GUID const& vessel_guid));
  MOCK_CONST_METHOD1(GetVessel, not_null<Vessel*>(GUID const& vessel_guid));

  not_null<std::unique_ptr<NavigationFrame>>
  NewBodyCentredNonRotatingNavigationFrame(
      Index reference_body_index) const override;
  not_null<std::unique_ptr<NavigationFrame>>
  NewBarycentricRotatingNavigationFrame(
      Index primary_index,
      Index secondary_index) const override;

  MOCK_CONST_METHOD2(
      FillBodyCentredNonRotatingNavigationFrame,
      void(Index reference_body_index,
           std::unique_ptr<NavigationFrame>* navigation_frame));
  MOCK_CONST_METHOD3(
      FillBarycentricRotatingNavigationFrame,
      void(Index primary_index,
           Index secondary_index,
           std::unique_ptr<NavigationFrame>* navigation_frame));

  // NOTE(phl): Needed because gMock 1.7.0 wants to copy the unique_ptr<>.
  void SetPlottingFrame(
    not_null<std::unique_ptr<NavigationFrame>> plotting_frame) override;

  MOCK_METHOD1(SetPlottingFrameConstRef,
               void(NavigationFrame const& plotting_frame));

  MOCK_CONST_METHOD0(GetPlottingFrame,
                     not_null<NavigationFrame const*>());

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

  MOCK_CONST_METHOD0(BarycentricToWorldSun,
                     OrthogonalMap<Barycentric, WorldSun>());

  MOCK_CONST_METHOD0(CurrentTime, Instant());

  MOCK_CONST_METHOD1(WriteToMessage,
                     void(not_null<serialization::Plugin*> message));
};

}  // namespace internal_plugin

using internal_plugin::MockPlugin;

}  // namespace ksp_plugin
}  // namespace principia
