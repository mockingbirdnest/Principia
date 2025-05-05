#include "ksp_plugin/planetarium.hpp"

#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "base/not_null.hpp"
#include "base/serialization.hpp"
#include "base/status_utilities.hpp"  // 🧙 For CHECK_OK.
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/perspective.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "geometry/signature.hpp"
#include "geometry/space.hpp"
#include "geometry/space_transformations.hpp"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "ksp_plugin/frames.hpp"
#include "numerics/quadrature.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/equipotential.hpp"
#include "physics/lagrange_equipotentials.hpp"
#include "physics/massive_body.hpp"
#include "physics/mock_continuous_trajectory.hpp"  // 🧙 For MockContinuousTrajectory.  // NOLINT
#include "physics/mock_ephemeris.hpp"  // 🧙 For MockEphemeris.
#include "physics/mock_rigid_reference_frame.hpp"  // 🧙 For MockRigidReferenceFrame.  // NOLINT
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "physics/rotating_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/serialization.hpp"
#include "testing_utilities/solar_system_factory.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace ksp_plugin {

using ::testing::_;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Le;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::SizeIs;
using namespace principia::base::_not_null;
using namespace principia::base::_serialization;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_perspective;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_sign;
using namespace principia::geometry::_signature;
using namespace principia::geometry::_space;
using namespace principia::geometry::_space_transformations;
using namespace principia::integrators::_embedded_explicit_runge_kutta_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_planetarium;
using namespace principia::numerics::_quadrature;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::physics::_continuous_trajectory;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_equipotential;
using namespace principia::physics::_lagrange_equipotentials;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_discrete_trajectory_factories;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_serialization;
using namespace principia::testing_utilities::_solar_system_factory;
using namespace principia::testing_utilities::_vanishes_before;

class PlanetariumTest : public ::testing::Test {
 protected:
  using LeftNavigation =
    Frame<struct LeftNavigationTag, Arbitrary, Handedness::Left>;

  PlanetariumTest()
      :  // The camera is located as {0, 20, 0} and is looking along -y.
        perspective_(
            RigidTransformation<Navigation, Camera>(
                Navigation::origin + Displacement<Navigation>(
                                         {0 * Metre, 20 * Metre, 0 * Metre}),
                Camera::origin,
                Rotation<LeftNavigation, Camera>(
                    Vector<double, LeftNavigation>({1, 0, 0}),
                    Vector<double, LeftNavigation>({0, 0, 1}),
                    Bivector<double, LeftNavigation>({0, -1, 0}))
                        .Forget<OrthogonalMap>() *
                    Signature<Navigation, LeftNavigation>(
                        Sign::Positive(),
                        Sign::Positive(),
                        DeduceSignReversingOrientation{})
                        .Forget<OrthogonalMap>())
                .Forget<Similarity>(),
            /*focal=*/5 * Metre),
        plotting_to_scaled_space_(
            [](Instant const&, Position<Navigation> const& plotted_point) {
              constexpr auto inverse_scale_factor = 1 / (6000 * Metre);
              return ScaledSpacePoint::FromCoordinates(
                  ((plotted_point - Navigation::origin) *
                   inverse_scale_factor).coordinates());
            }),
        // A body of radius 1 m located at the origin.
        body_(MassiveBody::Parameters(1 * Kilogram),
              RotatingBody<Barycentric>::Parameters(
                  /*mean_radius=*/1 * Metre,
                  /*reference_angle=*/0 * Radian,
                  /*reference_instant=*/t0_,
                  /*angular_frequency=*/10 * Radian / Second,
                  /*right_ascension_of_pole=*/0 * Radian,
                  /*declination_of_pole=*/π / 2 * Radian)),
        bodies_({&body_}),
        ephemeris_parameters_(
            SymmetricLinearMultistepIntegrator<
                QuinlanTremaine1990Order12,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
            /*step=*/10 * Minute),
        solar_system_(make_not_null_unique<SolarSystem<Barycentric>>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt",
            /*ignore_frame=*/true)),
        ephemeris_(solar_system_->MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            ephemeris_parameters_)),
        equipotential_parameters_(EmbeddedExplicitRungeKuttaIntegrator<
                                      DormandPrince1986RK547FC,
                                      Equipotential<Barycentric, World>::ODE>(),
                                  /*max_steps=*/1000,
                                  /*length_integration_tolerance=*/1 * Metre) {
    ON_CALL(plotting_frame_, t_min()).WillByDefault(Return(InfinitePast));
    ON_CALL(plotting_frame_, t_max()).WillByDefault(Return(InfiniteFuture));
    EXPECT_CALL(plotting_frame_, ToThisFrameAtTime(_))
        .WillRepeatedly(Return(RigidMotion<Barycentric, Navigation>(
            RigidTransformation<Barycentric, Navigation>::Identity(),
            Barycentric::nonrotating,
            Barycentric::unmoving)));
    EXPECT_CALL(mock_ephemeris_, bodies()).WillRepeatedly(ReturnRef(bodies_));
    EXPECT_CALL(mock_ephemeris_, trajectory(_))
        .WillRepeatedly(Return(&continuous_trajectory_));
    EXPECT_CALL(continuous_trajectory_, EvaluatePosition(_))
        .WillRepeatedly(Return(Barycentric::origin));
  }

  std::vector<std::unique_ptr<
      LagrangeEquipotentials<Barycentric, Navigation>::LinesBySpecificEnergy>>
  ComputeLagrangeEquipotentials() {
    auto const& earth = *solar_system_->massive_body(
        *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Earth));
    auto const& moon = *solar_system_->massive_body(
        *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Moon));

    std::vector<std::unique_ptr<
        LagrangeEquipotentials<Barycentric, Navigation>::LinesBySpecificEnergy>>
        result;

    LagrangeEquipotentials<Barycentric, Navigation> lagrange_equipotentials(
        ephemeris_.get());

    // Compute over 30 days.
    for (int i = 0; i < 30; ++i) {
      Instant const t = t0_ + i * Day;
      CHECK_OK(ephemeris_->Prolong(t));

      LagrangeEquipotentials<Barycentric, Navigation>::Parameters const
          equipotential_parameters{
              .primaries = {&earth}, .secondaries = {&moon}, .time = t};
      auto status_or_equipotentials =
          lagrange_equipotentials.ComputeLines(equipotential_parameters);
      CHECK_OK(status_or_equipotentials);
      result.emplace_back(
          std::make_unique<LagrangeEquipotentials<Barycentric, Navigation>::
                               LinesBySpecificEnergy>(
              std::move(status_or_equipotentials.value().lines)));
    }

    return result;
  }

  // It is convenient to represent the lines resulting from plotting as
  // discrete trajectories, even though they are not really physical, as it
  // makes it possible to evaluate them at any time.  The velocities are all
  // null.
  std::vector<DiscreteTrajectory<Navigation>> ComputePlottedLines(
      std::vector<std::unique_ptr<LagrangeEquipotentials<
          Barycentric,
          Navigation>::LinesBySpecificEnergy>> const& lagrange_equipotentials,
      Planetarium::Parameters const& planetarium_parameters,
      std::function<void(Planetarium const&,
                    DiscreteTrajectory<Navigation> const &)> const& plot) {
    auto const& earth = *solar_system_->massive_body(
        *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Earth));
    auto const& moon = *solar_system_->massive_body(
        *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Moon));

    auto const plotting_frame(
        RotatingPulsatingReferenceFrame<Barycentric, Navigation>(
            ephemeris_.get(), &earth, &moon));

    // The camera is located as {0, 0, 1} (one lunar distance above the
    // Earth-Moon barycentre) and is looking along -z.
    Perspective<Navigation, Camera> const perspective(
        RigidTransformation<Navigation, Camera>(
            Navigation::origin +
                Displacement<Navigation>({0 * Metre, 0 * Metre, 1 * Metre}),
            Camera::origin,
            Rotation<World, Camera>(Vector<double, World>({1, 0, 0}),
                                    Bivector<double, World>({0, -1, 0}),
                                    Vector<double, World>({0, 0, -1}))
                    .Forget<OrthogonalMap>() *
                Signature<Navigation, World>(
                    Sign::Positive(),
                    Signature<Navigation, World>::DeduceSign(),
                    Sign::Positive())
                    .Forget<OrthogonalMap>())
            .Forget<Similarity>(),
        /*focal=*/5 * Metre);

    std::vector<DiscreteTrajectory<Navigation>> plotted_trajectories;

    auto const plotting_to_scaled_space =
        [&plotted_trajectories](Instant const& plot_time,
                                Position<Navigation> const& plotted_point) {
          constexpr auto inverse_scale_factor = 1 / Metre;
          CHECK_OK(plotted_trajectories.back().Append(
              plot_time,
              DegreesOfFreedom<Navigation>(plotted_point,
                                           Velocity<Navigation>())));
          return ScaledSpacePoint::FromCoordinates(
              ((plotted_point - Navigation::origin) * inverse_scale_factor)
                  .coordinates());
        };

    Planetarium planetarium(planetarium_parameters,
                            perspective,
                            ephemeris_.get(),
                            &plotting_frame,
                            plotting_to_scaled_space);

    // Compute over 30 days.
    for (int i = 0; i < 30; ++i) {
      Instant const t = t0_ + i * Day;
      CHECK_OK(ephemeris_->Prolong(t));

      for (auto const& [_, lines] : *lagrange_equipotentials[i]) {
        for (auto const& line : lines) {
          plotted_trajectories.emplace_back();
          plot(planetarium, line);
        }
      }
    }

    return plotted_trajectories;
  }

  Length ComputePlottedLinesDistance(
      DiscreteTrajectory<Navigation> const& trajectory1,
      DiscreteTrajectory<Navigation> const& trajectory2) {
    CHECK_EQ(trajectory1.t_min(), trajectory2.t_min());
    CHECK_EQ(trajectory1.t_max(), trajectory2.t_max());
    return Sqrt(AutomaticClenshawCurtis(
                    [&trajectory1, &trajectory2](Instant const& t) {
                      return (trajectory1.EvaluatePosition(t) -
                              trajectory2.EvaluatePosition(t))
                          .Norm²();
                    },
                    trajectory1.t_min(),
                    trajectory2.t_max(),
                    /*max_relative_error=*/1.0e-6,
                    /*max_point=*/1000) /
                (trajectory1.t_max() - trajectory1.t_min()));
  }

  Instant const t0_;
  Perspective<Navigation, Camera> const perspective_;
  MockRigidReferenceFrame<Barycentric, Navigation> plotting_frame_;
  Planetarium::PlottingToScaledSpaceConversion plotting_to_scaled_space_;
  RotatingBody<Barycentric> const body_;
  std::vector<not_null<MassiveBody const*>> const bodies_;
  Ephemeris<Barycentric>::FixedStepParameters const ephemeris_parameters_;
  not_null<std::unique_ptr<SolarSystem<Barycentric>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
  Equipotential<Barycentric, World>::AdaptiveParameters const
      equipotential_parameters_;
  MockContinuousTrajectory<Barycentric> continuous_trajectory_;
  MockEphemeris<Barycentric> mock_ephemeris_;
};

TEST_F(PlanetariumTest, PlotMethod0) {
  DiscreteTrajectory<Barycentric> discrete_trajectory;
  AppendTrajectoryTimeline(/*from=*/NewCircularTrajectoryTimeline<Barycentric>(
                                        /*period=*/10 * Second,
                                        /*r=*/10 * Metre,
                                        /*Δt=*/1 * Second,
                                        /*t1=*/t0_,
                                        /*t2=*/t0_ + 11 * Second),
                           /*to=*/discrete_trajectory);

  // No dark area, infinite acuity, wide field of view.
  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1,
      /*angular_resolution=*/0 * Degree,
      /*field_of_view=*/90 * Degree);
  Planetarium planetarium(parameters,
                          perspective_,
                          &mock_ephemeris_,
                          &plotting_frame_,
                          plotting_to_scaled_space_);
  auto const rp2_lines =
      planetarium.PlotMethod0(discrete_trajectory,
                              discrete_trajectory.begin(),
                              discrete_trajectory.end(),
                              t0_ + 10 * Second,
                              /*reverse=*/false);

  // Because of the way the trajectory was constructed we have two lines which
  // meet in front of the camera and are separated by a hole behind the planet.
  EXPECT_THAT(rp2_lines, SizeIs(2));
  EXPECT_THAT(rp2_lines[0].front().x() - rp2_lines[1].back().x(),
              VanishesBefore(1 * Metre, 2));
  EXPECT_THAT(rp2_lines[0].back().x() - rp2_lines[1].front().x(),
              AlmostEquals(-10.0 / Sqrt(399.0) * Metre, 4, 94));

  for (auto const& rp2_line : rp2_lines) {
    for (auto const& rp2_point : rp2_line) {
      // The following limit is obtained by elementary geometry by noticing that
      // the trajectory is viewed from the camera under an angle of π / 6.
      EXPECT_THAT(rp2_point.x(),
                  AllOf(Ge(-5.0 / Sqrt(3.0) * Metre),
                        Le(5.0 / Sqrt(3.0) * Metre)));
      EXPECT_THAT(rp2_point.y(), VanishesBefore(1 * Metre, 0, 13));
    }
  }
}

TEST_F(PlanetariumTest, PlotMethod1) {
  // A quarter of a circular trajectory around the origin, with many small
  // segments.
  DiscreteTrajectory<Barycentric> discrete_trajectory;
  AppendTrajectoryTimeline(/*from=*/NewCircularTrajectoryTimeline<Barycentric>(
                                        /*period=*/100'000 * Second,
                                        /*r=*/10 * Metre,
                                        /*Δt=*/1 * Second,
                                        /*t1=*/t0_,
                                        /*t2=*/t0_ + 25'000 * Second),
                           /*to=*/discrete_trajectory);

  // No dark area, human visual acuity, wide field of view.
  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1,
      /*angular_resolution=*/0.4 * ArcMinute,
      /*field_of_view=*/90 * Degree);
  Planetarium planetarium(parameters,
                          perspective_,
                          &mock_ephemeris_,
                          &plotting_frame_,
                          plotting_to_scaled_space_);
  auto const rp2_lines =
      planetarium.PlotMethod1(discrete_trajectory,
                              discrete_trajectory.begin(),
                              discrete_trajectory.end(),
                              t0_ + 10 * Second,
                              /*reverse=*/false);

  EXPECT_THAT(rp2_lines, SizeIs(1));
  EXPECT_THAT(rp2_lines[0], SizeIs(4954));
  for (auto const& rp2_point : rp2_lines[0]) {
    EXPECT_THAT(rp2_point.x(),
                AllOf(Ge(0 * Metre),
                      Le(5.0 / Sqrt(3.0) * Metre)));
    EXPECT_THAT(rp2_point.y(), VanishesBefore(1 * Metre, 0, 14));
  }
}

TEST_F(PlanetariumTest, PlotMethod2) {
  // A quarter of a circular trajectory around the origin, with many small
  // segments.
  DiscreteTrajectory<Barycentric> discrete_trajectory;
  AppendTrajectoryTimeline(/*from=*/NewCircularTrajectoryTimeline<Barycentric>(
                                        /*period=*/100'000 * Second,
                                        /*r=*/10 * Metre,
                                        /*Δt=*/1 * Second,
                                        /*t1=*/t0_,
                                        /*t2=*/t0_ + 25'000 * Second),
                           /*to=*/discrete_trajectory);

  // No dark area, human visual acuity, wide field of view.
  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1,
      /*angular_resolution=*/0.4 * ArcMinute,
      /*field_of_view=*/90 * Degree);
  Planetarium planetarium(parameters,
                          perspective_,
                          &mock_ephemeris_,
                          &plotting_frame_,
                          plotting_to_scaled_space_);
  auto const rp2_lines =
      planetarium.PlotMethod2(discrete_trajectory,
                              discrete_trajectory.begin(),
                              discrete_trajectory.end(),
                              t0_ + 10 * Second,
                              /*reverse=*/false);

  EXPECT_THAT(rp2_lines, SizeIs(1));
  EXPECT_THAT(rp2_lines[0], SizeIs(43));
  for (auto const& rp2_point : rp2_lines[0]) {
    EXPECT_THAT(rp2_point.x(),
                AllOf(Ge(0 * Metre),
                      Le((5.0 / Sqrt(3.0)) * Metre)));
    EXPECT_THAT(rp2_point.y(), VanishesBefore(1 * Metre, 0, 14));
  }
}

#if !defined(_DEBUG)
TEST_F(PlanetariumTest, PlotMethod2_RealSolarSystem) {
  auto const discrete_trajectory =
      DiscreteTrajectory<Barycentric>::ReadFromMessage(
          ParseFromBytes<serialization::DiscreteTrajectory>(
              ReadFromBinaryFile(SOLUTION_DIR / "ksp_plugin_test" /
                                 "planetarium_trajectory.proto.bin")),
          /*tracked=*/{});

  auto ephemeris = Ephemeris<Barycentric>::ReadFromMessage(
      InfiniteFuture,
      ParseFromBytes<serialization::Ephemeris>(
          ReadFromBinaryFile(SOLUTION_DIR / "ksp_plugin_test" /
                             "planetarium_ephemeris.proto.bin")));

  auto plotting_frame = NavigationFrame::ReadFromMessage(
      ParseFromBytes<serialization::ReferenceFrame>(
          ReadFromBinaryFile(SOLUTION_DIR / "ksp_plugin_test" /
                             "planetarium_plotting_frame.proto.bin")),
      ephemeris.get());

  auto rigid_transformation =
      RigidTransformation<Navigation, Camera>::ReadFromMessage(
          ParseFromBytes<serialization::AffineMap>(
              ReadFromBinaryFile(SOLUTION_DIR / "ksp_plugin_test" /
                                 "planetarium_to_camera.proto.bin")));

  EXPECT_EQ(23423, discrete_trajectory.size());

  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1,
      /*angular_resolution=*/0.4 * ArcMinute,
      /*field_of_view=*/90 * Degree);
  Planetarium planetarium(
      parameters,
      Perspective<Navigation, Camera>(rigid_transformation.Forget<Similarity>(),
                                      /*focal=*/1 * Metre),
      ephemeris.get(),
      plotting_frame.get(),
      plotting_to_scaled_space_);
  auto const rp2_lines =
      planetarium.PlotMethod2(discrete_trajectory,
                              discrete_trajectory.begin(),
                              discrete_trajectory.end(),
                              discrete_trajectory.back().time,
                              /*reverse=*/false);

  EXPECT_EQ(2, rp2_lines.size());
  EXPECT_EQ(2, rp2_lines[0].size());
  EXPECT_EQ(9, rp2_lines[1].size());
}

TEST_F(PlanetariumTest, PlotMethod3_Equipotentials_AngularResolution) {
  // Human visual acuity.
  Angle reference_angular_resolution = 0.08 * ArcMinute;

  auto const lagrange_equipotentials = ComputeLagrangeEquipotentials();

  auto const plot_method3 =
      [](Planetarium const& planetarium,
         DiscreteTrajectory<Navigation> const& line) {
        planetarium.PlotMethod3(
            line,
            line.front().time,
            line.back().time,
            /*reverse=*/false,
            [](ScaledSpacePoint const&) {},
            /*max_points=*/std::numeric_limits<int>::max());
      };

  auto const plotted_trajectories_reference = ComputePlottedLines(
      lagrange_equipotentials,
      Planetarium::Parameters(
          /*sphere_radius_multiplier=*/1.0,
          reference_angular_resolution,
          /*field_of_view=*/90 * Degree),
      plot_method3);

  std::int64_t number_of_points_reference = 0;
  for (auto const& plotted_trajectory : plotted_trajectories_reference) {
    number_of_points_reference += plotted_trajectory.size();
  }
  EXPECT_EQ(100363, number_of_points_reference);

  // Compute plotted lines at progressively worse resolution and their distance
  // to the reference line.
  Angle trial_angular_resolution = reference_angular_resolution;
  for (std::int64_t i = 0; i < 8; ++i) {
    trial_angular_resolution *= 2;

    auto const plotted_trajectories =
        ComputePlottedLines(lagrange_equipotentials,
                            Planetarium::Parameters(
                                /*sphere_radius_multiplier=*/1.0,
                                trial_angular_resolution,
                                /*field_of_view=*/90 * Degree),
                            plot_method3);

    std::int64_t number_of_points = 0;
    for (auto const& plotted_trajectory : plotted_trajectories) {
      number_of_points += plotted_trajectory.size();
    }

    Length max_distance;
    for (std::int64_t i = 0; i < plotted_trajectories_reference.size(); ++i) {
      max_distance = std::max(
          max_distance,
          ComputePlottedLinesDistance(plotted_trajectories_reference[i],
                                      plotted_trajectories[i]));
    }

    switch (i) {
      case 0:
        // 0.8 arc min.
        EXPECT_EQ(71897, number_of_points);
        EXPECT_THAT(max_distance, IsNear(0.67_(1) * Centi(Metre)));
        break;
      case 1:
        // 1.6 arc min.
        EXPECT_EQ(51570, number_of_points);
        EXPECT_THAT(max_distance, IsNear(0.83_(1) * Centi(Metre)));
        break;
      case 2:
        // 3.2 arc min.
        EXPECT_EQ(36989, number_of_points);
        EXPECT_THAT(max_distance, IsNear(1.10_(1) * Centi(Metre)));
        break;
      case 3:
        // 6.4 arc min.
        EXPECT_EQ(26565, number_of_points);
        EXPECT_THAT(max_distance, IsNear(1.48_(1) * Centi(Metre)));
        break;
      case 4:
        // 12.8 arc min.
        EXPECT_EQ(19172, number_of_points);
        EXPECT_THAT(max_distance, IsNear(2.07_(1) * Centi(Metre)));
        break;
      case 5:
        // 25.6 arc min.
        EXPECT_EQ(13733, number_of_points);
        EXPECT_THAT(max_distance, IsNear(3.35_(1) * Centi(Metre)));
        break;
      case 6:
        // 51.2 arc min.
        EXPECT_EQ(9831, number_of_points);
        EXPECT_THAT(max_distance, IsNear(4.50_(1) * Centi(Metre)));
        break;
      case 7:
        // 102.4 arc min.
        EXPECT_EQ(7117, number_of_points);
        EXPECT_THAT(max_distance, IsNear(7.10_(1) * Centi(Metre)));
        break;
    }
  }
}
#endif

}  // namespace ksp_plugin
}  // namespace principia
