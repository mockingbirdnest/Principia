#include "physics/lagrange_equipotentials.hpp"

#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/body_centred_body_direction_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "physics/solar_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {

using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::Le;
using ::testing::Lt;
using ::testing::SizeIs;
using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::physics::_body_centred_body_direction_reference_frame;  // NOLINT
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_lagrange_equipotentials;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_solar_system_factory;

class LagrangeEquipotentialsTest : public ::testing::Test {
 protected:
  using Barycentric = Frame<struct BarycentricTag, Inertial>;
  using World = Frame<struct WorldTag, Arbitrary>;
  LagrangeEquipotentialsTest()
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
            ephemeris_parameters_)) {}

  Instant const t0_;
  Ephemeris<Barycentric>::FixedStepParameters const ephemeris_parameters_;
  not_null<std::unique_ptr<SolarSystem<Barycentric>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
};

#if !_DEBUG
TEST_F(LagrangeEquipotentialsTest,
       DISABLED_RotatingPulsating_GlobalOptimization) {
  Logger logger(TEMP_DIR / "equipotential_rp_global.wl",
                /*make_unique=*/false);
  std::int64_t const number_of_days = 502;
  auto const earth = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Earth));
  auto const moon = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Moon));
  auto const reference_frame(
      RotatingPulsatingReferenceFrame<Barycentric, World>(
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
      .periapsis_distance = 71'000 * Kilo(Metre),
      .apoapsis_distance = 0.65 * moon_elements.periapsis_distance.value(),
      .inclination = moon_elements.inclination,
      .longitude_of_ascending_node = moon_elements.longitude_of_ascending_node,
      .argument_of_periapsis = *moon_elements.argument_of_periapsis + Degree,
      .mean_anomaly = 0 * Degree};
  auto const earth_world_dof =
      reference_frame.ToThisFrameAtTimeSimilarly(t0_)(earth_dof);
  auto const moon_world_dof =
      reference_frame.ToThisFrameAtTimeSimilarly(t0_)(moon_dof);
  Position<World> const q_earth = earth_world_dof.position();
  Position<World> const q_moon = moon_world_dof.position();
  Position<World> const initial_earth_moon_l5 =
      Barycentre({q_earth, q_moon}, {1.0, 1.0}) +
      (q_earth - q_moon).Norm() *
          Vector<double, World>({0, quantities::Sqrt(3) / 2, 0});
  using MEO = Frame<struct MEOTag, Arbitrary>;
  BodyCentredBodyDirectionReferenceFrame<Barycentric, MEO> meo(
      ephemeris_.get(), moon, earth);
  // The initial states for four trajectories:
  // [0]: initially stationary in the rotating-pulsating frame near L3;
  // [1]: initially stationary in MEO at L5;
  // [2]: initially stationary in the rotating-pulsating frame at L5;
  // [3]: in an elliptic Earth orbit that reaches 65% of the way to the Moon.
  std::vector<DegreesOfFreedom<Barycentric>> const initial_states{
      reference_frame.FromThisFrameAtTimeSimilarly(t0_)(
          {q_earth + (q_earth - q_moon), World::unmoving}),
      meo.FromThisFrameAtTime(t0_)(
          {meo.ToThisFrameAtTime(t0_).rigid_transformation()(
               reference_frame.FromThisFrameAtTimeSimilarly(t0_).similarity()(
                   initial_earth_moon_l5)),
           MEO::unmoving}),
      reference_frame.FromThisFrameAtTimeSimilarly(t0_)(
          {initial_earth_moon_l5, World::unmoving}),
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

  LOG(ERROR) << "Flowing trajectories";
  CHECK_OK(
      ephemeris_->FlowWithFixedStep(t0_ + number_of_days * Day, *instance));
  LOG(ERROR) << "Flowed";

  Instant t = t0_;

  std::vector<std::vector<Position<World>>> trajectory_positions(
      trajectories.size());
  for (int j = 0; j < number_of_days; ++j) {
    LOG(ERROR) << "Day #" << j;
    t = t0_ + j * Day;
    CHECK_OK(ephemeris_->Prolong(t));
    for (int i = 0; i < trajectories.size(); ++i) {
      DegreesOfFreedom<World> const dof =
          reference_frame.ToThisFrameAtTimeSimilarly(t)(
              trajectories[i]->EvaluateDegreesOfFreedom(t));
      trajectory_positions[i].push_back(dof.position());
    }
    auto const equipotentials =
        LagrangeEquipotentials<Barycentric, World>(ephemeris_.get())
            .ComputeLines(
                {.primaries = {earth}, .secondaries = {moon}, .time = t});
    CHECK_OK(equipotentials.status());

    std::vector<SpecificEnergy> maxima;
    std::vector<Position<World>> arg_maximorum;
    for (auto const& [maximum, arg_maximi] : equipotentials->maxima) {
      EXPECT_THAT((arg_maximi - World::origin).Norm(),
                  AllOf(Gt(0.980 * Metre), Lt(1.016 * Metre)));
      maxima.push_back(maximum);
      arg_maximorum.push_back(arg_maximi);
    }
    EXPECT_THAT(maxima, SizeIs(AllOf(Ge(3), Le(6))));
    logger.Append("maxima", maxima, ExpressIn(Metre, Second));
    logger.Append("argMaximorum", arg_maximorum, ExpressIn(Metre));

    std::vector<SpecificEnergy> energies;
    std::vector<std::vector<std::vector<Position<World>>>> equipotentials_at_t;
    for (auto const& [energy, lines] : equipotentials->lines) {
      energies.push_back(energy);
      std::vector<std::vector<Position<World>>>& equipotentials_at_energy =
          equipotentials_at_t.emplace_back();
      for (auto const& line : lines) {
        std::vector<Position<World>>& equipotential =
            equipotentials_at_energy.emplace_back();
        for (auto const& [_, dof] : line) {
          EXPECT_THAT((dof.position() - World::origin).Norm(),
                      AllOf(Gt(0.736 * Metre), Lt(1.322 * Metre)));
          equipotential.push_back(dof.position());
        }
      }
    }
    logger.Append("energies", energies, ExpressIn(Metre, Second));
    logger.Append("equipotentialsEarthMoonGlobalOptimization",
                  equipotentials_at_t,
                  ExpressIn(Metre));
  }
  std::vector<std::vector<Position<World>>> world_trajectories;
  for (auto const& trajectory : trajectories) {
    world_trajectories.emplace_back();
    for (auto const& [t, dof] : *trajectory) {
      world_trajectories.back().push_back(
          reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
              dof.position()));
    }
  }
  logger.Set("trajectories", world_trajectories, ExpressIn(Metre));
  logger.Set("trajectoryPositions", trajectory_positions, ExpressIn(Metre));
}
#endif

}  // namespace physics
}  // namespace principia
