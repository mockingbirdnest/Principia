#include "physics/apsides.hpp"

#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.
#include "testing_utilities/string_log_sink.hpp"

namespace principia {
namespace physics {

using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::HasSubstr;
using ::testing::IsEmpty;
using ::testing::SizeIs;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::physics::_apsides;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_rotating_body;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_discrete_trajectory_factories;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_string_log_sink;

class ApsidesTest : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag, Inertial>;
};

#if !defined(_DEBUG)

TEST_F(ApsidesTest, ComputeApsidesDiscreteTrajectory) {
  Instant const t0;
  GravitationalParameter const μ = SolarGravitationalParameter;
  auto const b = new MassiveBody(μ);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<World>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b));
  initial_state.emplace_back(World::origin, World::unmoving);

  Ephemeris<World> ephemeris(
      std::move(bodies),
      initial_state,
      t0,
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Metre,
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<World>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<World>::NewtonianMotionEquation>(),
          10 * Minute));

  Displacement<World> r(
      {1 * AstronomicalUnit, 2 * AstronomicalUnit, 3 * AstronomicalUnit});
  Length const r_norm = r.Norm();
  Velocity<World> v({4 * Kilo(Metre) / Second,
                     5 * Kilo(Metre) / Second,
                     6 * Kilo(Metre) / Second});
  Speed const v_norm = v.Norm();

  Time const T = 2 * π * Sqrt(-(Pow<3>(r_norm) * Pow<2>(μ) /
                                Pow<3>(r_norm * Pow<2>(v_norm) - 2 * μ)));
  Length const a = -r_norm * μ / (r_norm * Pow<2>(v_norm) - 2 * μ);

  DiscreteTrajectory<World> trajectory;
  EXPECT_OK(trajectory.Append(t0,
                              DegreesOfFreedom<World>(World::origin + r, v)));

  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<World>::NoIntrinsicAcceleration,
      t0 + 10 * JulianYear,
      Ephemeris<World>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<World>::NewtonianMotionEquation>(),
          std::numeric_limits<std::int64_t>::max(),
          1e-3 * Metre,
          1e-3 * Metre / Second),
      Ephemeris<World>::unlimited_max_ephemeris_steps));

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(*ephemeris.trajectory(b),
                 trajectory,
                 trajectory.begin(),
                 trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_points=*/std::numeric_limits<int>::max(),
                 apoapsides,
                 periapsides);

  std::optional<Instant> previous_time;
  std::map<Instant, DegreesOfFreedom<World>> all_apsides;
  for (auto const& [time, degrees_of_freedom] : apoapsides) {
    all_apsides.emplace(time, degrees_of_freedom);
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(T, 118, 2824));
    }
    previous_time = time;
  }

  previous_time = std::nullopt;
  for (auto const& [time, degrees_of_freedom] : periapsides) {
    all_apsides.emplace(time, degrees_of_freedom);
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(T, 134, 257));
    }
    previous_time = time;
  }

  EXPECT_EQ(6, all_apsides.size());

  previous_time = std::nullopt;
  std::optional<Position<World>> previous_position;
  for (auto const& [time, degrees_of_freedom] : all_apsides) {
    Position<World> const position = degrees_of_freedom.position();
    if (previous_time) {
      EXPECT_THAT(time - *previous_time,
                  AlmostEquals(0.5 * T, 103, 5098));
      EXPECT_THAT((position - *previous_position).Norm(),
                  AlmostEquals(2.0 * a, 0, 176));
    }
    previous_time = time;
    previous_position = position;
  }
}

TEST_F(ApsidesTest, ComputeApsidesDiscreteTrajectory_Circular) {
  Instant const t0;
  Instant const t1 = t0 + 1e-15 * Second;
  Instant const t2 = t1 + 3 * Second;
  Time const Δt = 1.0 / 128.0 * Second;

  DiscreteTrajectory<World> reference_trajectory;
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewMotionlessTrajectoryTimeline(World::origin, Δt, t1, t2),
      reference_trajectory);
  AppendTrajectoryTimeline(
      NewCircularTrajectoryTimeline<World>(/*period=*/1 * Second,
                                           /*r=*/1 * Metre,
                                           Δt,
                                           t1,
                                           t2),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(reference_trajectory,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_points=*/10,
                 apoapsides,
                 periapsides);

  // This is a "suspicious" apsis, located at the beginning of a time interval,
  // because the circular trajectory leads to ill-conditioning.
  auto it = apoapsides.begin();
  EXPECT_THAT(it->time,
              AnyOf(AlmostEquals(t1 + 4 * Δt, 0),    // Windows, Ubuntu.
                    AlmostEquals(t1 + 2 * Δt, 0)));  // macOS.

  RotatingBody<World> const body(
      1 * Kilogram,
      RotatingBody<World>::Parameters(
          /*min_radius=*/1 * Metre,
          /*mean_radius=*/1 * Metre,
          /*max_radius=*/1 * Metre,
          /*reference_angle=*/0 * Radian,
          /*reference_instant=*/t0,
          /*angular_frequency=*/2 * π * Radian / Second,
          /*right_ascension_of_pole=*/0 * Radian,
          /*declination_of_pole=*/π / 2 * Radian));

  // The apsides do not oscillate in altitude because of the ill-conditioning,
  // so we give up.  This used to fail, see #3925.
  StringLogSink log_warning(google::WARNING);
  const auto intervals = ComputeCollisionIntervals(body,
                                                   reference_trajectory,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals, IsEmpty());
  EXPECT_THAT(log_warning.string(), HasSubstr("Anomalous apsides"));
}

#endif

TEST_F(ApsidesTest, ComputeFirstCollision) {
  Instant const t0;
  DiscreteTrajectory<World> reference_trajectory;
  DiscreteTrajectory<World> vessel_trajectory;

  // At `t0` the vessel is inside the celestial, so we expect the collision at a
  // negative time.
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({0 * Metre, -1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t0,
          /*t1=*/t0 - 10 * Second,
          /*t2=*/t0 + 10 * Second),
      reference_trajectory);

  // Note that the trajectory is short enough that we only have to go to degree
  // 64 for the interpolant.
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({1 * Metre, -1 * Metre, 0 * Metre}),
              Velocity<World>({0 * Metre / Second,
                               -1 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t0,
          /*t1=*/t0 - 10 * Second,
          /*t2=*/t0 + 1.9 * Second),
      vessel_trajectory);

  RotatingBody<World> const body(
      1 * Kilogram,
      RotatingBody<World>::Parameters(
          /*min_radius=*/1 * Metre,
          /*mean_radius=*/2 * Metre,
          /*max_radius=*/3 * Metre,
          /*reference_angle=*/0 * Radian,
          /*reference_instant=*/t0,
          /*angular_frequency=*/2 * π * Radian / Second,
          /*right_ascension_of_pole=*/0 * Radian,
          /*declination_of_pole=*/π / 2 * Radian));

  // The celestial is infinite in the z direction and has four lobes in the x-y
  // plane.  Think of a LEGO® axle.
  auto radius = [](Angle const& latitude, Angle const& longitude) {
    return (Cos(4 * longitude) + 2) * Metre;
  };

  // The computations below were verified with Mathematica to the given
  // accuracy.
  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(reference_trajectory,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_points=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, IsEmpty());
  EXPECT_THAT(periapsides, SizeIs(1));
  EXPECT_THAT(periapsides.begin()->time, AlmostEquals(t0 + 0.5 * Second, 0));

  const auto intervals = ComputeCollisionIntervals(body,
                                                   reference_trajectory,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals,
              ElementsAre(IntervalMatches(
                  AlmostEquals(t0 + (1.0 - Sqrt(17.0)) / 2.0 * Second, 1),
                  AlmostEquals(t0 + 1 * Second, 0))));

  auto const maybe_collision =
      ComputeFirstCollision(body,
                            reference_trajectory,
                            vessel_trajectory,
                            intervals[0],
                            /*max_error=*/2e-4 * Metre,
                            radius);
  auto const& collision = maybe_collision.value();

  EXPECT_THAT(collision.time - t0,
              IsNear(-1.43862_(1) * Second));
  EXPECT_THAT(collision.degrees_of_freedom.position() - World::origin,
              Componentwise(1 * Metre,
                            IsNear(0.43862_(1) * Metre),
                            0 * Metre));
  EXPECT_THAT(
      collision.degrees_of_freedom.velocity(),
      AlmostEquals(Velocity<World>({0 * Metre / Second,
                                    -1 * Metre / Second,
                                    0 * Metre / Second}), 0));
}

#if !defined(_DEBUG)

TEST_F(ApsidesTest, ComputeNodes) {
  Instant const t0;
  GravitationalParameter const μ = SolarGravitationalParameter;
  auto const b = new MassiveBody(μ);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<World>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b));
  initial_state.emplace_back(World::origin, World::unmoving);

  Ephemeris<World> ephemeris(
      std::move(bodies),
      initial_state,
      t0,
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Metre,
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<World>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<World>::NewtonianMotionEquation>(),
          10 * Minute));

  KeplerianElements<World> elements;
  elements.eccentricity = 0.25;
  elements.semimajor_axis = 1 * AstronomicalUnit;
  elements.inclination = 10 * Degree;
  elements.longitude_of_ascending_node = 42 * Degree;
  elements.argument_of_periapsis = 100 * Degree;
  elements.mean_anomaly = 0 * Degree;
  KeplerOrbit<World> const orbit{
      *ephemeris.bodies()[0], MasslessBody{}, elements, t0};
  elements = orbit.elements_at_epoch();

  DiscreteTrajectory<World> trajectory;
  EXPECT_OK(trajectory.Append(t0, initial_state[0] + orbit.StateVectors(t0)));

  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<World>::NoIntrinsicAcceleration,
      t0 + 10 * JulianYear,
      Ephemeris<World>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<World>::NewtonianMotionEquation>(),
          std::numeric_limits<std::int64_t>::max(),
          1e-3 * Metre,
          1e-3 * Metre / Second),
      Ephemeris<World>::unlimited_max_ephemeris_steps));

  Vector<double, World> const north({0, 0, 1});

  DiscreteTrajectory<World> ascending_nodes;
  DiscreteTrajectory<World> descending_nodes;
  EXPECT_OK(ComputeNodes(trajectory,
                         trajectory.begin(),
                         trajectory.end(),
                         /*t_max=*/InfiniteFuture,
                         north,
                         /*max_points=*/std::numeric_limits<int>::max(),
                         ascending_nodes,
                         descending_nodes));

  std::optional<Instant> previous_time;
  for (auto const& [time, degrees_of_freedom] : ascending_nodes) {
    EXPECT_THAT((degrees_of_freedom.position() - World::origin)
                    .coordinates()
                    .ToSpherical()
                    .longitude,
                AlmostEquals(elements.longitude_of_ascending_node, 0, 104));
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(*elements.period, 0, 20));
    }
    previous_time = time;
  }

  previous_time = std::nullopt;
  for (auto const& [time, degrees_of_freedom] : descending_nodes) {
    EXPECT_THAT(
        (degrees_of_freedom.position() - World::origin)
                .coordinates()
                .ToSpherical()
                .longitude,
        AlmostEquals(elements.longitude_of_ascending_node - π * Radian, 0, 29));
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(*elements.period, 0, 29));
    }
    previous_time = time;
  }

  EXPECT_THAT(ascending_nodes, SizeIs(10));
  EXPECT_THAT(descending_nodes, SizeIs(10));

  DiscreteTrajectory<World> south_ascending_nodes;
  DiscreteTrajectory<World> south_descending_nodes;
  Vector<double, World> const mostly_south({1, 1, -1});
  EXPECT_OK(ComputeNodes(trajectory,
                         trajectory.begin(),
                         trajectory.end(),
                         /*t_max=*/InfiniteFuture,
                         mostly_south,
                         /*max_points=*/std::numeric_limits<int>::max(),
                         south_ascending_nodes,
                         south_descending_nodes));
  EXPECT_THAT(south_ascending_nodes, SizeIs(10));
  EXPECT_THAT(south_descending_nodes, SizeIs(10));

  for (auto south_ascending_it  = south_ascending_nodes.begin(),
            ascending_it        = ascending_nodes.begin(),
            south_descending_it = south_descending_nodes.begin(),
            descending_it       = descending_nodes.begin();
       south_ascending_it != south_ascending_nodes.end();
       ++south_ascending_it,
       ++ascending_it,
       ++south_descending_it,
       ++descending_it) {
    EXPECT_THAT(south_ascending_it->degrees_of_freedom,
                Eq(descending_it->degrees_of_freedom));
    EXPECT_THAT(south_ascending_it->time, Eq(descending_it->time));
    EXPECT_THAT(south_descending_it->degrees_of_freedom,
                Eq(ascending_it->degrees_of_freedom));
    EXPECT_THAT(south_descending_it->time, Eq(ascending_it->time));
  }
}

#endif

// A dedicated fixture for `ComputeCollisionIntervals` because we have many
// tests for that function.
class ApsidesTest_ComputeCollisionIntervals : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag, Inertial>;

  ApsidesTest_ComputeCollisionIntervals()
      :  // Only `max_radius` matters for the body.
        body_(1 * Kilogram,
              RotatingBody<World>::Parameters(
                  /*min_radius=*/2 * Metre,
                  /*mean_radius=*/2 * Metre,
                  /*max_radius=*/2 * Metre,
                  /*reference_angle=*/0 * Radian,
                  /*reference_instant=*/t0_,
                  /*angular_frequency=*/2 * π * Radian / Second,
                  /*right_ascension_of_pole=*/0 * Radian,
                  /*declination_of_pole=*/π / 2 * Radian)) {
    // The celestial is motionless at origin.
    AppendTrajectoryTimeline(NewMotionlessTrajectoryTimeline(World::origin,
                                                             /*Δt=*/1 * Second,
                                                             t1_,
                                                             t2_),
                             body_trajectory_);
  }

  Instant const t0_;
  Instant const t1_ = t0_;
  Instant const t2_ = t0_ + 20 * Second;
  RotatingBody<World> const body_;
  DiscreteTrajectory<World> body_trajectory_;
};

// A linear trajectory that intersects the body.  There is one periapsis below
// `max_radius` and two sentinel apoapsides at the extremities of the
// trajectory.
TEST_F(ApsidesTest_ComputeCollisionIntervals, OnePeriapsisBelowMaxRadius) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({-4 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t2_),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, IsEmpty());
  EXPECT_THAT(periapsides, SizeIs(1));

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals,
              ElementsAre(IntervalMatches(
                  AlmostEquals(t0_ + (4 - Sqrt(3)) * Second, 0),
                  AlmostEquals(t0_ + (4 + Sqrt(3)) * Second, 0))));
}

// A linear trajectory that does not intersect the body.  There is one periapsis
// above `max_radius`.
TEST_F(ApsidesTest_ComputeCollisionIntervals, OnePeriapsisAboveMaxRadius) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({-4 * Metre, 3 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t2_),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, IsEmpty());
  EXPECT_THAT(periapsides, SizeIs(1));

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals, IsEmpty());
}

// A linear trajectory without a periapsis or an apoapsis.
TEST_F(ApsidesTest_ComputeCollisionIntervals, NoPeriapsis) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({-30 * Metre, 3 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t2_),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, IsEmpty());
  EXPECT_THAT(periapsides, IsEmpty());

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals, IsEmpty());
}

// Two perpendicular linear trajectory segments resulting in an apoapsis
// below `max_radius` and two periapsides.
TEST_F(ApsidesTest_ComputeCollisionIntervals, OneApoapsisBelowMaxRadius) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({-4 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t0_ + 4.9 * Second),
      vessel_trajectory);
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({1 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({0 * Metre / Second,
                               -1 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t0_ + 5 * Second,
          t2_),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, SizeIs(1));
  EXPECT_THAT(periapsides, SizeIs(2));

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals,
              ElementsAre(IntervalMatches(
                  AlmostEquals(t0_ + (4 - Sqrt(3)) * Second, 0),
                  AlmostEquals(t0_ + (6 + Sqrt(3)) * Second, 1))));
}

// Two linear trajectory segments at 45° resulting in an apoapsis above
// `max_radius` and two periapsides.
TEST_F(ApsidesTest_ComputeCollisionIntervals, OneApoapsisAboveMaxRadius) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({-4 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t0_ + 6.9 * Second),
      vessel_trajectory);
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({3 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({-1 * Metre / Second,
                               -1 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t0_ + 7 * Second,
          t2_),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, SizeIs(1));
  EXPECT_THAT(periapsides, SizeIs(2));

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals,
              ElementsAre(IntervalMatches(
                              AlmostEquals(t0_ + (4 - Sqrt(3)) * Second, 0),
                              AlmostEquals(t0_ + (4 + Sqrt(3)) * Second, 1)),
                          IntervalMatches(AlmostEquals(t0_ + 8 * Second, 0),
                                          AlmostEquals(t0_ + 10 * Second, 0))));
}

// Two linear trajectory segments at 45° resulting in an apoapsis above
// `max_radius` followed by a periapsis below `max_radius`.  The initial point
// acts as a periapsis below `max_radius`.
TEST_F(ApsidesTest_ComputeCollisionIntervals,
       InitialPeriapsisOneApoapsisOnePeriapsis) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({1 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t0_ + 1.9 * Second),
      vessel_trajectory);
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({3 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({-1 * Metre / Second,
                               -1 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t0_ + 2 * Second,
          t2_),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, SizeIs(1));
  EXPECT_THAT(periapsides, SizeIs(1));

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals,
              ElementsAre(IntervalMatches(
                              AlmostEquals(t0_, 0),
                              AlmostEquals(t0_ + (Sqrt(3) - 1) * Second, 2)),
                          IntervalMatches(AlmostEquals(t0_ + 3 * Second, 0),
                                          AlmostEquals(t0_ + 5 * Second, 0))));
}

// Two linear trajectory segments at 45° resulting in a periapsis below
// `max_radius` followed by an apoapsis above `max_radius`.  The final point
// acts as a periapsis below `max_radius`.
TEST_F(ApsidesTest_ComputeCollisionIntervals,
       OnePeriapsisOneApoapsisFinalPeriapsis) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({-3 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t0_ + 5.9 * Second),
      vessel_trajectory);
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({3 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({-0.9 * Metre / Second,
                               -0.9 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t0_ + 6 * Second,
          t0_ + 8.9 * Second),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, SizeIs(1));
  EXPECT_THAT(periapsides, SizeIs(1));

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(
      intervals,
      ElementsAre(
          IntervalMatches(AlmostEquals(t0_ + (3 - Sqrt(3)) * Second, 1),
                          AlmostEquals(t0_ + (3 + Sqrt(3)) * Second, 1)),
          IntervalMatches(AlmostEquals(t0_ + (6 + 1 / 0.9) * Second, 1),
                          AlmostEquals(t0_ + 8 * Second, 0))));
}

// Two linear trajectory segments at 45° resulting in a periapsis below
// `max_radius` followed by an apoapsis above `max_radius`.
TEST_F(ApsidesTest_ComputeCollisionIntervals, OnePeriapsisOneApoapsis) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({-4 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t0_ + 6.9 * Second),
      vessel_trajectory);
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({3 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({-1 * Metre / Second,
                               -1 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t0_ + 7 * Second,
          t0_ + 8.9 * Second),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, SizeIs(1));
  EXPECT_THAT(periapsides, SizeIs(1));

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals,
              ElementsAre(IntervalMatches(
                  AlmostEquals(t0_ + (4 - Sqrt(3)) * Second, 0),
                  AlmostEquals(t0_ + (4 + Sqrt(3)) * Second, 1))));
}

// A linear trajectory with a periapsis below `max_radius`.  The initial point
// acts as an apoapsis below `max_radius`.
TEST_F(ApsidesTest_ComputeCollisionIntervals, InitialApoapsisOnePeriapsis) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({-1 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t2_),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, IsEmpty());
  EXPECT_THAT(periapsides, SizeIs(1));

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals,
              ElementsAre(IntervalMatches(
                              AlmostEquals(t0_, 0),
                              AlmostEquals(t0_ + (1 + Sqrt(3)) * Second, 0))));
}

// A linear trajectory with a periapsis below `max_radius`.  The final point
// acts as an apoapsis below `max_radius`.
TEST_F(ApsidesTest_ComputeCollisionIntervals, OnePeriapsisFinalApoapsis) {
  DiscreteTrajectory<World> vessel_trajectory;
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({-18 * Metre, 1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t1_,
          t2_),
      vessel_trajectory);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(body_trajectory_,
                 vessel_trajectory,
                 vessel_trajectory.begin(),
                 vessel_trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_point=*/10,
                 apoapsides,
                 periapsides);
  EXPECT_THAT(apoapsides, IsEmpty());
  EXPECT_THAT(periapsides, SizeIs(1));

  const auto intervals = ComputeCollisionIntervals(body_,
                                                   body_trajectory_,
                                                   vessel_trajectory,
                                                   apoapsides,
                                                   periapsides);
  EXPECT_THAT(intervals,
              ElementsAre(IntervalMatches(
                              AlmostEquals(t0_ + (18 - Sqrt(3)) * Second, 0),
                              AlmostEquals(t0_ + 19 * Second, 0))));
}

}  // namespace physics
}  // namespace principia
