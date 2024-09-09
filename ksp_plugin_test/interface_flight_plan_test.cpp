#include "ksp_plugin/interface.hpp"

#include <memory>
#include <string>
#include <utility>

#include "base/not_null.hpp"
#include "geometry/identity.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/permutation.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin/renderer.hpp"
#include "ksp_plugin/vessel.hpp"
#include "ksp_plugin_test/mock_flight_plan.hpp"  // 🧙 For MockFlightPlan.
#include "ksp_plugin_test/mock_manœuvre.hpp"  // 🧙 For MockManœuvre.
#include "ksp_plugin_test/mock_plugin.hpp"  // 🧙 For MockPlugin.
#include "ksp_plugin_test/mock_renderer.hpp"  // 🧙 For MockRenderer.
#include "ksp_plugin_test/mock_vessel.hpp"  // 🧙 For MockVessel.
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/mock_continuous_trajectory.hpp"  // 🧙 For MockContinuousTrajectory.  // NOLINT
#include "physics/mock_ephemeris.hpp"  // 🧙 For MockEphemeris.
#include "physics/mock_rigid_reference_frame.hpp"  // 🧙 For MockRigidReferenceFrame.  // NOLINT
#include "physics/rigid_motion.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "quantities/constants.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace interface {

using ::testing::AllOf;
using ::testing::AnyNumber;
using ::testing::ByMove;
using ::testing::DoAll;
using ::testing::Invoke;
using ::testing::Property;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::SetArgReferee;
using ::testing::StrictMock;
using ::testing::_;
using namespace principia::base::_not_null;
using namespace principia::geometry::_identity;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_permutation;
using namespace principia::geometry::_rotation;
using namespace principia::integrators::_embedded_explicit_generalized_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_manœuvre;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin::_renderer;
using namespace principia::ksp_plugin::_vessel;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_continuous_trajectory;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::quantities::_constants;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

namespace {

char const vessel_guid[] = "123-456";
Index const celestial_index = 1;

MATCHER_P(HasThrust, thrust, "") {
  return arg.thrust == thrust;
}

MATCHER_P(HasSpecificImpulse, specific_impulse, "") {
  return arg.specific_impulse == specific_impulse;
}

MATCHER_P(HasInitialTime, initial_time, "") {
  return arg.timing.initial_time && *arg.timing.initial_time == initial_time;
}

MATCHER_P(HasΔv, Δv, "") {
  return arg.intensity.Δv && *arg.intensity.Δv == Δv;
}

MATCHER(IsOk,
        std::string(negation ? "is not" : "is") + " ok") {
  return arg.error == 0;
}

}  // namespace

class InterfaceFlightPlanTest : public ::testing::Test {
 protected:
  InterfaceFlightPlanTest()
      : adaptive_step_parameters_(
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
            /*max_steps=*/111,
            /*length_integration_tolerance=*/222 * Metre,
            /*speed_integration_tolerance=*/333 * Metre / Second),
        generalized_adaptive_step_parameters_(
            EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                Fine1987RKNG34,
                Ephemeris<Barycentric>::GeneralizedNewtonianMotionEquation>(),
            /*max_steps=*/111,
            /*length_integration_tolerance=*/222 * Metre,
            /*speed_integration_tolerance=*/333 * Metre / Second),
        identity_(Rotation<Barycentric, AliceSun>::Identity()),
        plugin_(make_not_null_unique<StrictMock<MockPlugin>>()),
        const_plugin_(plugin_.get()) {}

  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters_;
  Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
      generalized_adaptive_step_parameters_;
  Rotation<Barycentric, AliceSun> identity_;
  std::unique_ptr<MockManœuvre<Barycentric, Navigation>> navigation_manœuvre_;
  MockRenderer renderer_;
  StrictMock<MockFlightPlan> flight_plan_;
  not_null<std::unique_ptr<StrictMock<MockPlugin>>> const plugin_;
  StrictMock<MockPlugin> const* const const_plugin_;
  Instant const t0_;
};

TEST_F(InterfaceFlightPlanTest, FlightPlan) {
  Burn interface_burn = {
      /*thrust_in_kilonewtons=*/1,
      /*specific_impulse_in_seconds_g0=*/2,
      /*frame=*/{/*extension=*/6000, /*centre=*/celestial_index},
      /*initial_time=*/3,
      /*delta_v=*/{4, 5, 6},
      /*is_inertially_fixed=*/true};
  StrictMock<MockVessel> vessel;

  EXPECT_CALL(*plugin_, HasVessel(vessel_guid))
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*plugin_, GetVessel(vessel_guid))
      .WillRepeatedly(Return(&vessel));
  EXPECT_CALL(vessel, has_flight_plan())
      .WillRepeatedly(Return(true));
  EXPECT_CALL(vessel, flight_plan())
      .WillRepeatedly(ReturnRef(flight_plan_));
  EXPECT_CALL(*plugin_, ExtendPredictionForFlightPlan(vessel_guid))
      .Times(AnyNumber());

  EXPECT_TRUE(principia__FlightPlanExists(plugin_.get(), vessel_guid));

  EXPECT_CALL(*plugin_, CreateFlightPlan(vessel_guid,
                                         Instant() + 30 * Second,
                                         100 * Tonne));
  principia__FlightPlanCreate(plugin_.get(),
                              vessel_guid,
                              /*final_time=*/30,
                              /*mass_in_tonnes=*/100);

  EXPECT_CALL(flight_plan_, SetDesiredFinalTime(Instant() + 60 * Second))
      .WillOnce(Return(absl::OkStatus()));
  EXPECT_THAT(
      *principia__FlightPlanSetDesiredFinalTime(plugin_.get(), vessel_guid, 60),
      IsOk());

  EXPECT_CALL(flight_plan_, initial_time())
      .WillOnce(Return(Instant() + 3 * Second));
  EXPECT_EQ(3, principia__FlightPlanGetInitialTime(plugin_.get(), vessel_guid));

  EXPECT_CALL(flight_plan_, desired_final_time())
      .WillOnce(Return(Instant() + 4 * Second));
  EXPECT_EQ(4, principia__FlightPlanGetDesiredFinalTime(plugin_.get(),
                                                        vessel_guid));

  EXPECT_CALL(
      flight_plan_,
      SetAdaptiveStepParameters(
          AllOf(
              Property(
                  &Ephemeris<Barycentric>::AdaptiveStepParameters::max_steps,
                  11),
              Property(&Ephemeris<Barycentric>::AdaptiveStepParameters::
                            length_integration_tolerance,
                        22 * Metre),
              Property(&Ephemeris<Barycentric>::AdaptiveStepParameters::
                            speed_integration_tolerance<>,
                        33 * Metre / Second)),
          AllOf(
              Property(&Ephemeris<Barycentric>::
                            GeneralizedAdaptiveStepParameters::max_steps,
                        11),
              Property(
                  &Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters::
                      length_integration_tolerance,
                  22 * Metre),
              Property(
                  &Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters::
                      speed_integration_tolerance<>,
                  33 * Metre / Second))))
      .WillOnce(Return(absl::OkStatus()));
  EXPECT_THAT(*principia__FlightPlanSetAdaptiveStepParameters(
                  plugin_.get(),
                  vessel_guid,
                  {/*integrator_kind=*/1,
                   /*generalized_integrator_kind=*/2,
                   /*max_step=*/11,
                   /*length_integration_tolerance=*/22,
                   /*speed_integration_tolerance=*/33}),
              IsOk());

  EXPECT_CALL(flight_plan_, adaptive_step_parameters())
      .WillOnce(ReturnRef(adaptive_step_parameters_));
  EXPECT_CALL(flight_plan_, generalized_adaptive_step_parameters())
      .WillOnce(ReturnRef(generalized_adaptive_step_parameters_));
  FlightPlanAdaptiveStepParameters expected_adaptive_step_parameters = {
      /*integrator_kind=*/1,
      /*generalized_integrator_kind=*/2,
      /*max_step=*/111,
      /*length_integration_tolerance=*/222,
      /*speed_integration_tolerance=*/333};
  EXPECT_EQ(expected_adaptive_step_parameters,
            principia__FlightPlanGetAdaptiveStepParameters(
                plugin_.get(), vessel_guid));

  EXPECT_CALL(*plugin_,
              NewBodyCentredNonRotatingNavigationFrame(celestial_index))
      .WillOnce(Return(
          ByMove(std::make_unique<StrictMock<
                     MockRigidReferenceFrame<Barycentric, Navigation>>>())));
  EXPECT_CALL(
      flight_plan_,
      Insert(AllOf(HasThrust(1 * Kilo(Newton)),
                   HasSpecificImpulse(2 * Second * StandardGravity),
                   HasInitialTime(Instant() + 3 * Second),
                   HasΔv(Velocity<Frenet<Navigation>>({4 * (Metre / Second),
                                                       5 * (Metre / Second),
                                                       6 * (Metre / Second)}))),
             0))
      .WillOnce(Return(absl::OkStatus()));
  EXPECT_THAT(*principia__FlightPlanInsert(
                  plugin_.get(), vessel_guid, interface_burn, 0),
              IsOk());

  EXPECT_CALL(flight_plan_, number_of_manœuvres())
      .WillOnce(Return(4));
  EXPECT_EQ(4, principia__FlightPlanNumberOfManoeuvres(plugin_.get(),
                                                       vessel_guid));

  auto const plotting_frame =
      make_not_null_unique<MockRigidReferenceFrame<Barycentric, Navigation>>();

  MockEphemeris<Barycentric> ephemeris;
  MassiveBody const centre(MassiveBody::Parameters("centre", 1 * Kilogram));
  MockContinuousTrajectory<Barycentric> centre_trajectory;
  EXPECT_CALL(ephemeris, trajectory(check_not_null(&centre)))
      .WillOnce(Return(&centre_trajectory));
  // Cannot use a mock here since we use `dynamic_cast` to find the type of the
  // actual frame.
  BodyCentredNonRotatingReferenceFrame<Barycentric, Navigation> const* const
      navigation_manœuvre_frame =
          new BodyCentredNonRotatingReferenceFrame<Barycentric, Navigation>(
            &ephemeris,
            &centre);

  MockManœuvre<Barycentric, Navigation>::Intensity intensity;
  intensity.direction = Vector<double, Frenet<Navigation>>({1, 1, 1});
  intensity.duration = 7 * Second;
  MockManœuvre<Barycentric, Navigation>::Timing timing;
  timing.initial_time = Instant();
  MockManœuvre<Barycentric, Navigation>::Burn const burn{
      intensity,
      timing,
      10 * Kilo(Newton),
      30 * Second * StandardGravity,
      std::unique_ptr<RigidReferenceFrame<Barycentric, Navigation> const>(
          navigation_manœuvre_frame),
      /*is_inertially_fixed=*/true};
  navigation_manœuvre_ =
      std::make_unique<MockManœuvre<Barycentric, Navigation>>(20 * Tonne, burn);
  auto const barycentric_to_plotting = RigidMotion<Barycentric, Navigation>(
      RigidTransformation<Barycentric, Navigation>(
          Barycentric::origin,
          Navigation::origin,
          OrthogonalMap<Barycentric, Navigation>::Identity()),
      Barycentric::nonrotating,
      Barycentric::unmoving);
  EXPECT_CALL(*plugin_, renderer()).WillRepeatedly(ReturnRef(renderer_));
  EXPECT_CALL(*const_plugin_, renderer()).WillRepeatedly(ReturnRef(renderer_));
  EXPECT_CALL(*plugin_, PlanetariumRotation())
      .WillRepeatedly(ReturnRef(identity_));
  EXPECT_CALL(flight_plan_, GetManœuvre(3))
      .WillOnce(ReturnRef(*navigation_manœuvre_));
  EXPECT_CALL(*plugin_, CelestialIndexOfBody(Ref(centre)))
      .WillOnce(Return(celestial_index));
  auto const navigation_manoeuvre =
      std::unique_ptr<NavigationManoeuvre>(
          principia__FlightPlanGetManoeuvre(plugin_.get(),
                                            vessel_guid,
                                            3));

  EXPECT_EQ(10, navigation_manoeuvre->burn.thrust_in_kilonewtons);
  EXPECT_EQ(6000, navigation_manoeuvre->burn.frame.extension);
  EXPECT_EQ(celestial_index, navigation_manoeuvre->burn.frame.centre_index);
  EXPECT_EQ(20, navigation_manoeuvre->initial_mass_in_tonnes);
  EXPECT_THAT(navigation_manoeuvre->burn.specific_impulse_in_seconds_g0,
              AlmostEquals(30, 1));

  EXPECT_CALL(flight_plan_, GetManœuvre(3))
      .WillOnce(ReturnRef(*navigation_manœuvre_));

  EXPECT_CALL(renderer_, BarycentricToWorldSun(_))
      .WillOnce(Return(
          Permutation<Barycentric, WorldSun>(
              Permutation<Barycentric, WorldSun>::CoordinatePermutation::YXZ)
              .Forget<OrthogonalMap>()));
  EXPECT_CALL(*navigation_manœuvre_, FrenetFrame())
      .WillOnce(
          Return(OrthogonalMap<Frenet<Navigation>, Barycentric>::Identity()));
  EXPECT_CALL(*plugin_, CurrentTime())
      .WillRepeatedly(Return(Instant() - 4 * Second));
  EXPECT_CALL(renderer_, GetPlottingFrame())
      .WillRepeatedly(Return(plotting_frame.get()));
  EXPECT_CALL(*plotting_frame, ToThisFrameAtTime(Instant()))
      .WillOnce(Return(barycentric_to_plotting));
  EXPECT_CALL(*plotting_frame, FromThisFrameAtTime(Instant() - 4 * Second))
      .WillOnce(Return(barycentric_to_plotting.Inverse()));
  principia__FlightPlanGetManoeuvreFrenetTrihedron(plugin_.get(),
                                                   vessel_guid,
                                                   3);

  EXPECT_CALL(flight_plan_, number_of_segments())
      .WillOnce(Return(12));
  EXPECT_EQ(12, principia__FlightPlanNumberOfSegments(plugin_.get(),
                                                      vessel_guid));

  DiscreteTrajectory<World> rendered_trajectory;
  EXPECT_OK(rendered_trajectory.Append(
      t0_, DegreesOfFreedom<World>(World::origin, World::unmoving)));
  EXPECT_OK(rendered_trajectory.Append(
      t0_ + 1 * Second,
      DegreesOfFreedom<World>(
          World::origin +
              Displacement<World>({0 * Metre, 1 * Metre, 2 * Metre}),
          World::unmoving)));
  EXPECT_OK(rendered_trajectory.Append(
      t0_ + 2 * Second,
      DegreesOfFreedom<World>(
          World::origin +
              Displacement<World>({0 * Metre, 2 * Metre, 4 * Metre}),
          World::unmoving)));
  DiscreteTrajectory<Barycentric> segment;
  DegreesOfFreedom<Barycentric> immobile_origin{Barycentric::origin,
                                                Barycentric::unmoving};
  EXPECT_OK(segment.Append(t0_, immobile_origin));
  EXPECT_OK(segment.Append(t0_ + 1 * Second, immobile_origin));
  EXPECT_OK(segment.Append(t0_ + 2 * Second, immobile_origin));
  EXPECT_CALL(flight_plan_, GetSegment(3))
      .WillOnce(Return(segment.segments().begin()));
  EXPECT_CALL(renderer_, RenderBarycentricTrajectoryInWorld(_, _, _, _, _))
      .WillOnce(Return(ByMove(std::move(rendered_trajectory))));
  auto* const iterator =
      principia__FlightPlanRenderedSegment(plugin_.get(),
                                           vessel_guid,
                                           {0, 1, 2},
                                           3);
  EXPECT_EQ(XYZ({0, 0, 0}),
            principia__IteratorGetDiscreteTrajectoryXYZ(iterator));
  principia__IteratorIncrement(iterator);
  EXPECT_EQ(XYZ({0, 1, 2}),
            principia__IteratorGetDiscreteTrajectoryXYZ(iterator));
  principia__IteratorIncrement(iterator);
  EXPECT_EQ(XYZ({0, 2, 4}),
            principia__IteratorGetDiscreteTrajectoryXYZ(iterator));

  interface_burn.thrust_in_kilonewtons = 10;
  EXPECT_CALL(*plugin_,
              NewBodyCentredNonRotatingNavigationFrame(celestial_index))
      .WillOnce(Return(
          ByMove(std::make_unique<StrictMock<
                     MockRigidReferenceFrame<Barycentric, Navigation>>>())));
  auto const manœuvre = NavigationManœuvre(/*initial_mass=*/1 * Kilogram, burn);
  EXPECT_CALL(
      flight_plan_,
      Replace(
          AllOf(HasThrust(10 * Kilo(Newton)),
                HasSpecificImpulse(2 * Second * StandardGravity),
                HasInitialTime(Instant() + 3 * Second),
                HasΔv(Velocity<Frenet<Navigation>>({4 * (Metre / Second),
                                                    5 * (Metre / Second),
                                                    6 * (Metre / Second)}))),
          42))
      .WillOnce(Return(absl::OkStatus()));
  EXPECT_THAT(*principia__FlightPlanReplace(
                  plugin_.get(), vessel_guid, interface_burn, 42),
              IsOk());

  EXPECT_CALL(flight_plan_, Remove(0));
  principia__FlightPlanRemove(plugin_.get(), vessel_guid, 0);

  EXPECT_CALL(vessel, DeleteFlightPlan());
  principia__FlightPlanDelete(plugin_.get(), vessel_guid);
}

}  // namespace interface
}  // namespace principia
