#include "physics/equipotential.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/plane.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/logger.hpp"
#include "numerics/global_optimization.hpp"
#include "numerics/root_finders.hpp"
#include "physics/body_centred_body_direction_dynamic_frame.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {

using integrators::EmbeddedExplicitRungeKuttaIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::DormandPrince1986RK547FC;
using integrators::methods::QuinlanTremaine1990Order12;
using integrators::methods::Quinlan1999Order8A;
using testing_utilities::SolarSystemFactory;
using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_named_quantities;
using namespace principia::geometry::_plane;
using namespace principia::geometry::_rotation;
using namespace principia::numerics::_global_optimization;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

class EquipotentialTest : public ::testing::Test {
 protected:
  using Barycentric = Frame<struct BarycentricTag, Inertial>;
  using World = Frame<struct WorldTag, Arbitrary>;

  EquipotentialTest()
      : ephemeris_parameters_(
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
        equipotential_parameters_(
            EmbeddedExplicitRungeKuttaIntegrator<
                DormandPrince1986RK547FC,
                Equipotential<Barycentric, World>::ODE>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Metre) {}

  Position<World> ComputePositionInWorld(
      Instant const& t,
      DynamicFrame<Barycentric, World> const& dynamic_frame,
      SolarSystemFactory::Index const body) {
    auto const to_this_frame = dynamic_frame.ToThisFrameAtTime(t);
    return to_this_frame.rigid_transformation()(
        solar_system_->trajectory(*ephemeris_, SolarSystemFactory::name(body))
            .EvaluatePosition(t));
  }

  std::array<Position<World>, 2> ComputeLagrangePoints(
      SolarSystemFactory::Index const body1,
      SolarSystemFactory::Index const body2,
      Instant const& t,
      DynamicFrame<Barycentric, World> const& dynamic_frame,
      Plane<World> const& plane) {
    auto const body1_position =
        ComputePositionInWorld(t, dynamic_frame, body1);
    auto const body2_position =
        ComputePositionInWorld(t, dynamic_frame, body2);
    auto const body2_body1 = body1_position - body2_position;

    auto const binormal = plane.UnitBinormals().front();
    Rotation<World, World> const rot_l4(-60 * Degree, binormal);
    auto const body2_l4 = rot_l4(body2_body1);
    auto const l4 = body2_l4 + body2_position;
    Rotation<World, World> const rot_l5(60 * Degree, binormal);
    auto const body2_l5 = rot_l5(body2_body1);
    auto const l5 = body2_l5 + body2_position;

    return {l4, l5};
  }

  // Logs to Mathematica the equipotential line for the given |body| in the
  // specified |dynamic_frame|.
  void LogEquipotentialLine(
      mathematica::Logger& logger,
      Plane<World> const& plane,
      Instant const& t,
      DynamicFrame<Barycentric, World> const& dynamic_frame,
      SolarSystemFactory::Index const body,
      std::string_view const suffix = "") {
    Equipotential<Barycentric, World> const equipotential(
        equipotential_parameters_, &dynamic_frame);
    std::string const name = SolarSystemFactory::name(body);

    CHECK_OK(ephemeris_->Prolong(t));
    auto const line =
        equipotential.ComputeLine(
            plane, t, ComputePositionInWorld(t0_, dynamic_frame, body));
    std::vector<Position<World>> positions;
    std::vector<double> βs;
    for (auto const& state : line) {
      auto const& [position, β] = state;
      positions.push_back(position);
      βs.push_back(β);
    }
    logger.Set(absl::StrCat("equipotential", name, suffix),
               positions,
               mathematica::ExpressIn(Metre));
    logger.Set(absl::StrCat("beta", name, suffix), βs);
  }

  // Logs to Mathematica a family of equipotential lines determined by a
  // parameter.  There must exist an overload of |ComputeLine| with a
  // |LineParameter| as its third argument.
  template<typename LineParameter>
  void LogFamilyOfEquipotentialLines(
      mathematica::Logger& logger,
      DynamicFrame<Barycentric, World> const& dynamic_frame,
      int const number_of_days,
      std::string_view const suffix,
      std::function<std::vector<LineParameter>(
          Position<World> const& l4,
          Position<World> const& l5)> const& get_line_parameters) {
    Equipotential<Barycentric, World> const equipotential(
        equipotential_parameters_, &dynamic_frame);
    auto const plane =
        Plane<World>::OrthogonalTo(Vector<double, World>({0, 0, 1}));

    std::vector<std::vector<std::vector<Position<World>>>> all_positions;
    std::vector<std::vector<std::vector<double>>> all_βs;
    for (int j = 0; j < number_of_days; ++j) {
      Instant const t = t0_ + j * Day;
      CHECK_OK(ephemeris_->Prolong(t));
      all_positions.emplace_back();
      all_βs.emplace_back();

      auto const& [l4, l5] = ComputeLagrangePoints(SolarSystemFactory::Earth,
                                                    SolarSystemFactory::Moon,
                                                    t,
                                                    dynamic_frame,
                                                    plane);

      for (auto const& line_parameter : get_line_parameters(l4, l5)) {
        auto const line =
            equipotential.ComputeLine(plane, t, line_parameter);
        all_positions.back().emplace_back();
        all_βs.back().emplace_back();
        for (auto const& state : line) {
          auto const& [position, β] = state;
          all_positions.back().back().push_back(position);
          all_βs.back().back().push_back(β);
        }
      }
    }
    logger.Set(absl::StrCat("equipotentialsEarthMoon", suffix),
               all_positions,
               mathematica::ExpressIn(Metre));
    logger.Set(absl::StrCat("betasEarthMoon", suffix), all_βs);
  }

  Instant const t0_;
  Ephemeris<Barycentric>::FixedStepParameters const ephemeris_parameters_;
  not_null<std::unique_ptr<SolarSystem<Barycentric>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
  Equipotential<Barycentric, World>::AdaptiveParameters const
      equipotential_parameters_;
};

#if !_DEBUG
TEST_F(EquipotentialTest, BodyCentredNonRotating) {
  mathematica::Logger logger(TEMP_DIR / "equipotential_bcnr.wl",
                             /*make_unique=*/false);
  auto const dynamic_frame(
      BodyCentredNonRotatingDynamicFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Sun))));
  Equipotential<Barycentric, World> const equipotential(
      equipotential_parameters_, &dynamic_frame);

  auto const plane =
      Plane<World>::OrthogonalTo(Vector<double, World>({2, 3, -5}));

  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       dynamic_frame,
                       SolarSystemFactory::Mercury);
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       dynamic_frame,
                       SolarSystemFactory::Earth);
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       dynamic_frame,
                       SolarSystemFactory::Jupiter, "Close");
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 100 * Day,
                       dynamic_frame,
                       SolarSystemFactory::Jupiter, "Far");
}

TEST_F(EquipotentialTest, BodyCentredBodyDirection_EquidistantPoints) {
  mathematica::Logger logger(TEMP_DIR / "equipotential_bcbd_distances.wl",
                             /*make_unique=*/false);
  auto const dynamic_frame(
      BodyCentredBodyDirectionDynamicFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Earth)),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Moon))));

  LogFamilyOfEquipotentialLines<Position<World>>(
      logger,
      dynamic_frame,
      /*number_of_days=*/30,
      /*suffix=*/"Distances",
      [](Position<World> const& l4, Position<World> const& l5) {
        std::vector<Position<World>> positions;
        for (int i = 0; i <= 10; ++i) {
          positions.push_back(Barycentre(
              std::pair{l4, l5}, std::pair{i / 10.0, (10.0 - i) / 10.0}));
        }
        return positions;
      });
}

TEST_F(EquipotentialTest, BodyCentredBodyDirection_EquidistantEnergies) {
  mathematica::Logger logger(TEMP_DIR / "equipotential_bcbd_energies.wl",
                             /*make_unique=*/false);
  auto const dynamic_frame(
      BodyCentredBodyDirectionDynamicFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Earth)),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Moon))));

  LogFamilyOfEquipotentialLines<DegreesOfFreedom<World>>(
      logger,
      dynamic_frame,
      /*number_of_days=*/30,
      /*suffix=*/"Energies",
      [](Position<World> const& l4, Position<World> const& l5) {
        auto const midpoint = Barycentre(std::pair{l4, l5}, std::pair{1, 1});
        std::vector<DegreesOfFreedom<World>> degrees_of_freedom;
        for (int i = 0; i <= 10; ++i) {
          degrees_of_freedom.push_back(
              DegreesOfFreedom<World>(midpoint,
                                      Velocity<World>({i * 100 * Metre / Second,
                                                       0 * Metre / Second,
                                                       0 * Metre / Second})));
        }
        for (int i = 0; i < 4; ++i) {
          degrees_of_freedom.push_back(DegreesOfFreedom<World>(
              midpoint,
              Velocity<World>({(1100 + i * 10) * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})));
        }
        return degrees_of_freedom;
      });
}

TEST_F(EquipotentialTest, BodyCentredBodyDirection_GlobalOptimization) {
  mathematica::Logger logger(TEMP_DIR / "equipotential_bcbd_global.wl",
                             /*make_unique=*/false);
  std::int64_t const number_of_days = 100;
  auto const earth = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Earth));
  auto const moon = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Moon));
  auto const dynamic_frame(
      BodyCentredBodyDirectionDynamicFrame<Barycentric, World>(
          ephemeris_.get(), moon, earth));
  CHECK_OK(ephemeris_->Prolong(t0_ + number_of_days * Day));

  DegreesOfFreedom<Barycentric> const earth_dof =
      ephemeris_->trajectory(earth)->EvaluateDegreesOfFreedom(t0_);
  DegreesOfFreedom<Barycentric> const moon_dof =
      ephemeris_->trajectory(moon)->EvaluateDegreesOfFreedom(t0_);
  KeplerOrbit<Barycentric> const moon_orbit(
      *earth,
      *moon,
      moon_dof - earth_dof,
      t0_);
  KeplerianElements<Barycentric> const moon_elements =
      moon_orbit.elements_at_epoch();

  KeplerianElements<Barycentric> const elements{
      .periapsis_distance = 7000 * Kilo(Metre),
      .apoapsis_distance = 0.8 * *moon_elements.periapsis_distance,
      .inclination = moon_elements.inclination,
      .longitude_of_ascending_node = moon_elements.longitude_of_ascending_node,
      .argument_of_periapsis = *moon_elements.argument_of_periapsis + Degree,
      .mean_anomaly = 0 * Degree};
  auto const earth_world_dof = dynamic_frame.ToThisFrameAtTime(t0_)(earth_dof);
  auto const moon_world_dof = dynamic_frame.ToThisFrameAtTime(t0_)(moon_dof);
  Position<World> const q_earth = earth_world_dof.position();
  Position<World> const q_moon = moon_world_dof.position();
  Velocity<World> const v_earth = earth_world_dof.velocity();
  Velocity<World> const v_moon = moon_world_dof.velocity();
  Position<World> const initial_earth_moon_l5 =
      Barycentre(std::pair(q_earth, q_moon), std::pair(1, 1)) +
      (q_earth - q_moon).Norm() *
          Vector<double, World>({0, quantities::Sqrt(3) / 2, 0});
  // The initial states for three trajectories:
  // [0]: initially stationary in |dynamic_frame| at L5;
  // [1]: initially stationary in the rotating-pulsating frame at L5;
  // [2]: in an elliptic Earth orbit that reaches 80% of the way to the Moon.
  std::vector<DegreesOfFreedom<Barycentric>> const initial_states{
      dynamic_frame.FromThisFrameAtTime(t0_)(
          {q_moon + 2 * (q_earth - q_moon), World::unmoving}),
      dynamic_frame.FromThisFrameAtTime(t0_)(
          {initial_earth_moon_l5, World::unmoving}),
      dynamic_frame.FromThisFrameAtTime(t0_)(
          {initial_earth_moon_l5,
           Barycentre(std::pair(v_earth, v_moon), std::pair(1, 1)) +
               InnerProduct(q_earth - q_moon, v_earth - v_moon) /
                   (q_earth - q_moon).Norm() *
                   Vector<double, World>({0, quantities::Sqrt(3) / 2, 0})}),
      earth_dof +
          KeplerOrbit<Barycentric>(*earth, MasslessBody{}, elements, t0_)
              .StateVectors(t0_)};

  std::vector<not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>>>
      trajectories;
  std::vector<not_null<DiscreteTrajectory<Barycentric>*>> instance_trajectories;
  for (auto const& s : initial_states) {
    trajectories.push_back(
        make_not_null_unique<DiscreteTrajectory<Barycentric>>());
    instance_trajectories.push_back(trajectories.back().get());
    CHECK_OK(trajectories.back()->Append(t0_, s));
    trajectories.back()->segments().front().SetDownsampling(
        {.max_dense_intervals = 10'000, .tolerance = 10 * Metre});
  }
  auto const instance = ephemeris_->NewInstance(
      instance_trajectories,
      Ephemeris<Barycentric>::NoIntrinsicAccelerations,
      Ephemeris<Barycentric>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              Quinlan1999Order8A,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*step=*/10 * Second));

  CHECK_OK(
      ephemeris_->FlowWithFixedStep(t0_ + number_of_days * Day, *instance));

  Instant t = t0_;
  auto const potential = [&dynamic_frame,
                          &t](Position<World> const& position) {
    return dynamic_frame.GeometricPotential(t, position);
  };
  auto const acceleration = [&dynamic_frame,
                              &t](Position<World> const& position) {
    auto const acceleration =
        dynamic_frame.GeometricAcceleration(t, {position, Velocity<World>{}});
    // Note the sign.
    return -Vector<Acceleration, World>({acceleration.coordinates()[0],
                                         acceleration.coordinates()[1],
                                         Acceleration{}});
  };
  const MultiLevelSingleLinkage<SpecificEnergy, Position<World>, 2>::Box box = {
      .centre = World::origin,
      .vertices = {Displacement<World>({1'000'000 * Kilo(Metre),
                                        0 * Metre,
                                        0 * Metre}),
                   Displacement<World>({0 * Metre,
                                        1'000'000 * Kilo(Metre),
                                        0 * Metre})}};

  MultiLevelSingleLinkage<SpecificEnergy, Position<World>, 2> optimizer(
      box, potential, acceleration);
  Equipotential<Barycentric, World> const equipotential(
      equipotential_parameters_, &dynamic_frame);
  auto const plane =
      Plane<World>::OrthogonalTo(Vector<double, World>({0, 0, 1}));

  std::vector<std::vector<std::vector<std::vector<Position<World>>>>>
      all_positions;
  std::vector<std::vector<std::vector<std::vector<double>>>> all_βs;
  std::vector<std::vector<Position<World>>> trajectory_positions(
      trajectories.size());
  std::vector<SpecificEnergy> energies;
  for (int j = 0; j < number_of_days; ++j) {
    LOG(ERROR) << "Day #" << j;
    t = t0_ + j * Day;
    CHECK_OK(ephemeris_->Prolong(t));
    all_positions.emplace_back();
    all_βs.emplace_back();

    auto const arg_maximorum = optimizer.FindGlobalMaxima(
        /*points_per_round=*/1000,
        /*number_of_rounds=*/std::nullopt,
        /*local_search_tolerance=*/10'000 * Kilo(Metre));
    logger.Append("argMaximorum", arg_maximorum, mathematica::ExpressIn(Metre));
    std::vector<SpecificEnergy> maxima;
    SpecificEnergy maximum_maximorum = -Infinity<SpecificEnergy>;

    Position<World> const earth_position =
        dynamic_frame.ToThisFrameAtTime(t).rigid_transformation()(
            ephemeris_->trajectory(earth)->EvaluatePosition(t));
    Position<World> const moon_position =
        dynamic_frame.ToThisFrameAtTime(t).rigid_transformation()(
            ephemeris_->trajectory(moon)->EvaluatePosition(t));
    for (auto const& arg_maximum : arg_maximorum) {
      maxima.push_back(potential(arg_maximum));
      maximum_maximorum = std::max(maximum_maximorum, maxima.back());
    }
    logger.Append("maxima", maxima, mathematica::ExpressIn(Metre, Second));

    double const arg_approx_l1 = numerics::Brent(
        [&](double x) {
          return potential(Barycentre(std::pair(moon_position, earth_position),
                                      std::pair(x, 1 - x)));
        },
        0.0,
        1.0,
        std::greater<>{});
    SpecificEnergy const approx_l1_energy =
        potential(Barycentre(std::pair(moon_position, earth_position),
                             std::pair(arg_approx_l1, 1 - arg_approx_l1)));

    for (int i = 0; i < trajectories.size(); ++i) {
      DegreesOfFreedom<World> const dof = dynamic_frame.ToThisFrameAtTime(t)(
          trajectories[i]->EvaluateDegreesOfFreedom(t));
      trajectory_positions[i].push_back(dof.position());
    }
    SpecificEnergy const ΔV = maximum_maximorum - approx_l1_energy;
    for (int i = 1; i <= 8; ++i) {
      all_positions.back().emplace_back();
      all_βs.back().emplace_back();
      SpecificEnergy const energy = maximum_maximorum - i * (1.0 / 7.0 * ΔV);
      auto const& lines = equipotential.ComputeLines(
          plane,
          t,
          arg_maximorum,
          {{moon_position, moon->min_radius()},
           {earth_position, earth->min_radius()}},
          [](Position<World> q) {
            return World::origin +
                   Normalize(q - World::origin) * 2'000'000 * Kilo(Metre);
          },
          energy);
      for (auto const& line : lines) {
        all_positions.back().back().emplace_back();
        all_βs.back().back().emplace_back();
        for (auto const& state : line) {
          auto const& [position, β] = state;
          all_positions.back().back().back().push_back(position);
          all_βs.back().back().back().push_back(β);
        }
      }
    }
  }
  std::vector<std::vector<Position<World>>> world_trajectories;
  for (auto const& trajectory : trajectories) {
    world_trajectories.emplace_back();
    for (auto const& [t, dof] : *trajectory) {
      world_trajectories.back().push_back(
          dynamic_frame.ToThisFrameAtTime(t).rigid_transformation()(
              dof.position()));
    }
  }
  std::vector<std::vector<Vector<double, World>>> pulsating_trajectories;
  for (auto const& trajectory : trajectories) {
    pulsating_trajectories.emplace_back();
    for (auto const& [t, dof] : *trajectory) {
      pulsating_trajectories.back().push_back(
          (dynamic_frame.ToThisFrameAtTime(t).rigid_transformation()(
               dof.position()) -
           World::origin) /
          (ephemeris_->trajectory(earth)->EvaluatePosition(t) -
           ephemeris_->trajectory(moon)->EvaluatePosition(t))
              .Norm());
    }
  }
  logger.Set("trajectories", world_trajectories, mathematica::ExpressIn(Metre));
  logger.Set("pulsatingTrajectories", pulsating_trajectories);
  logger.Set("trajectoryPositions",
             trajectory_positions,
             mathematica::ExpressIn(Metre));
  logger.Set("energies",
             energies,
             mathematica::ExpressIn(Metre, Second));
  logger.Set("equipotentialsEarthMoonGlobalOptimization",
             all_positions,
             mathematica::ExpressIn(Metre));
  logger.Set("betasEarthMoonGlobalOptimization", all_βs);
}

#endif

}  // namespace physics
}  // namespace principia
