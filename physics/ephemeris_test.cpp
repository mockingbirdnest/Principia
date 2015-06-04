#include "physics/ephemeris.hpp"

#include <map>
#include <vector>

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system.hpp"

namespace principia {

using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::Abs;
using quantities::ArcTan;
using quantities::Area;
using quantities::Pow;
using quantities::Sqrt;
using si::Kilogram;
using si::Metre;
using si::Milli;
using si::Minute;
using si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::ICRFJ2000Ecliptic;
using testing_utilities::kSolarSystemBarycentre;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystem;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;

namespace physics {

class EphemerisTest : public testing::Test {
 protected:
  using EarthMoonOrbitPlane = Frame<serialization::Frame::TestTag,
                                    serialization::Frame::TEST, true>;

  void SetUpEarthMoonSystem(
      not_null<std::vector<not_null<std::unique_ptr<MassiveBody const>>>*> const
          bodies,
      not_null<std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>>*> const
          initial_state,
      not_null<Position<EarthMoonOrbitPlane>*> const centre_of_mass,
      not_null<Time*> const period) {
    auto earth = std::make_unique<MassiveBody const>(6E24 * Kilogram);
    auto moon = std::make_unique<MassiveBody const>(7E22 * Kilogram);

    // The Earth-Moon system, roughly, with a circular orbit with velocities
    // in the centre-of-mass frame.
    Position<EarthMoonOrbitPlane> const q1(
        Vector<Length, EarthMoonOrbitPlane>({0 * Metre, 0 * Metre, 0 * Metre}));
    Position<EarthMoonOrbitPlane> const q2(
        Vector<Length, EarthMoonOrbitPlane>({0 * Metre,
                                             4E8 * Metre,
                                             0 * Metre}));
    Length const semi_major_axis = (q1 - q2).Norm();
    *period = 2 * π * Sqrt(Pow<3>(semi_major_axis) /
                               (earth->gravitational_parameter() +
                                moon->gravitational_parameter()));
    *centre_of_mass =
        geometry::Barycentre<Vector<Length, EarthMoonOrbitPlane>, Mass>(
            {q1, q2}, {earth->mass(), moon->mass()});
    Velocity<EarthMoonOrbitPlane> const v1(
        {-2 * π * (q1 - *centre_of_mass).Norm() / *period,
         0 * SIUnit<Speed>(),
         0 * SIUnit<Speed>()});
    Velocity<EarthMoonOrbitPlane> const v2(
        {2 * π * (q2 - *centre_of_mass).Norm() / *period,
         0 * SIUnit<Speed>(),
         0 * SIUnit<Speed>()});

    bodies->push_back(std::move(earth));
    bodies->push_back(std::move(moon));
    initial_state->push_back(DegreesOfFreedom<EarthMoonOrbitPlane>(q1, v1));
    initial_state->push_back(DegreesOfFreedom<EarthMoonOrbitPlane>(q2, v2));
  }

  Instant t0_;
};

// The canonical Earth-Moon system, tuned to produce circular orbits.
TEST_F(EphemerisTest, EarthMoon) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          period / 100,
          0.1 * Milli(Metre),
          5 * Milli(Metre));

  ephemeris.Prolong(t0_ + period);

  ContinuousTrajectory<EarthMoonOrbitPlane> const& earth_trajectory =
      ephemeris.trajectory(earth);
  ContinuousTrajectory<EarthMoonOrbitPlane> const& moon_trajectory =
      ephemeris.trajectory(moon);

  ContinuousTrajectory<EarthMoonOrbitPlane>::Hint hint;
  std::vector<Displacement<EarthMoonOrbitPlane>> earth_positions;
  for (int i = 0; i <= 100; ++i) {
    earth_positions.push_back(
      earth_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
          centre_of_mass);
  }
  EXPECT_THAT(earth_positions.size(), Eq(101));
  EXPECT_THAT(Abs(earth_positions[25].coordinates().y), Lt(6E-4 * Metre));
  EXPECT_THAT(Abs(earth_positions[50].coordinates().x), Lt(6E-3 * Metre));
  EXPECT_THAT(Abs(earth_positions[75].coordinates().y), Lt(2E-2 * Metre));
  EXPECT_THAT(Abs(earth_positions[100].coordinates().x), Lt(3E-2 * Metre));

  std::vector<Displacement<EarthMoonOrbitPlane>> moon_positions;
  for (int i = 0; i <= 100; ++i) {
    moon_positions.push_back(
      moon_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
          centre_of_mass);
  }
  EXPECT_THAT(moon_positions.size(), Eq(101));
  EXPECT_THAT(Abs(moon_positions[25].coordinates().y), Lt(5E-2 * Metre));
  EXPECT_THAT(Abs(moon_positions[50].coordinates().x), Lt(6E-1 * Metre));
  EXPECT_THAT(Abs(moon_positions[75].coordinates().y), Lt(2 * Metre));
  EXPECT_THAT(Abs(moon_positions[100].coordinates().x), Lt(2 * Metre));
}

// The Moon alone.  It moves in straight line.
TEST_F(EphemerisTest, Moon) {
  Position<EarthMoonOrbitPlane> const reference_position;
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  bodies.erase(bodies.begin());
  initial_state.erase(initial_state.begin());

  MassiveBody const* const moon = bodies[0].get();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          period / 100,
          0.1 * Milli(Metre),
          5 * Milli(Metre));

  ephemeris.Prolong(t0_ + period);

  ContinuousTrajectory<EarthMoonOrbitPlane> const& moon_trajectory =
      ephemeris.trajectory(moon);

  ContinuousTrajectory<EarthMoonOrbitPlane>::Hint hint;
  DegreesOfFreedom<EarthMoonOrbitPlane> const moon_degrees_of_freedom =
      moon_trajectory.EvaluateDegreesOfFreedom(t0_ + period, &hint);
  Length const q = (moon_degrees_of_freedom.position() -
                    reference_position).coordinates().y;
  Speed const v = moon_degrees_of_freedom.velocity().coordinates().x;
  std::vector<Displacement<EarthMoonOrbitPlane>> moon_positions;
  for (int i = 0; i <= 100; ++i) {
    moon_positions.push_back(
        moon_trajectory.EvaluatePosition(t0_ + i * period / 100, &hint) -
            reference_position);
  }

  EXPECT_THAT(moon_positions.size(), Eq(101));
  EXPECT_THAT(moon_positions[25].coordinates().x,
              AlmostEquals(0.25 * period * v, 13));
  EXPECT_THAT(moon_positions[25].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[50].coordinates().x,
              AlmostEquals(0.50 * period * v, 14));
  EXPECT_THAT(moon_positions[50].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[75].coordinates().x,
              AlmostEquals(0.75 * period * v, 19));
  EXPECT_THAT(moon_positions[75].coordinates().y, Eq(q));
  EXPECT_THAT(moon_positions[100].coordinates().x,
              AlmostEquals(1.00 * period * v, 13));
  EXPECT_THAT(moon_positions[100].coordinates().y, Eq(q));
}

// The Earth and a massless probe 1 billion meters away, with the same velocity,
// and an acceleration which exactly compensates gravitational attraction.  Both
// bodies move in straight lines.
TEST_F(EphemerisTest, EarthProbe) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody const* const earth = bodies[0].get();
  MassiveBody const* const moon = bodies[1].get();

  Ephemeris<EarthMoonOrbitPlane>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0_,
          McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
          period / 100,
          0.1 * Milli(Metre),
          5 * Milli(Metre));

  Trajectory<EarthMoonOrbitPlane> trajectory(earth);
  ephemeris.Flow(&trajectory,
                 0 * Metre,
                 0 * Metre / Second,
                 DormandElMikkawyPrince1986RKN434FM<
                     Position<EarthMoonOrbitPlane>>(),
                 t0_ + 1 * Second);

  //Position<EarthMoonOrbitPlane> const reference_position =
  //    Position<EarthMoonOrbitPlane>();
  //Length const distance = 1E9 * SIUnit<Length>();
  //trajectory3_->Append(trajectory1_->last().time(),
  //                     {trajectory1_->last().degrees_of_freedom().position() +
  //                          Vector<Length, EarthMoonOrbitPlane>(
  //                              {0 * SIUnit<Length>(),
  //                               distance,
  //                               0 * SIUnit<Length>()}),
  //                      trajectory1_->last().degrees_of_freedom().velocity()});
  //trajectory3_->set_intrinsic_acceleration(
  //    [this, distance](Instant const& t) {
  //  return Vector<Acceleration, EarthMoonOrbitPlane>(
  //      {0 * SIUnit<Acceleration>(),
  //       body1_.gravitational_parameter() / (distance * distance),
  //       0 * SIUnit<Acceleration>()});});

  //system_->Integrate(*integrator_,
  //                   trajectory1_->last().time() + period_,
  //                   period_ / 100,
  //                   1,      // sampling_period
  //                   false,  // tmax_is_exact
  //                   {trajectory1_.get(), trajectory3_.get()});

  //Length const q1 = (trajectory1_->last().degrees_of_freedom().position() -
  //                   reference_position).coordinates().y;
  //Speed const v1 =
  //    trajectory1_->last().degrees_of_freedom().velocity().coordinates().x;
  //std::vector<Vector<Length, EarthMoonOrbitPlane>> const positions1 =
  //    ValuesOf(trajectory1_->Positions(), reference_position);
  //LOG(INFO) << ToMathematicaString(positions1);
  //EXPECT_THAT(positions1.size(), Eq(101));
  //EXPECT_THAT(positions1[25].coordinates().x,
  //            AlmostEquals(0.25 * period_ * v1, 1));
  //EXPECT_THAT(positions1[25].coordinates().y, Eq(q1));
  //EXPECT_THAT(positions1[50].coordinates().x,
  //            AlmostEquals(0.50 * period_ * v1, 1));
  //EXPECT_THAT(positions1[50].coordinates().y, Eq(q1));
  //EXPECT_THAT(positions1[75].coordinates().x,
  //            AlmostEquals(0.75 * period_ * v1, 1, 2));
  //EXPECT_THAT(positions1[75].coordinates().y, Eq(q1));
  //EXPECT_THAT(positions1[100].coordinates().x,
  //            AlmostEquals(1.00 * period_ * v1, 1));
  //EXPECT_THAT(positions1[100].coordinates().y, Eq(q1));

  //Length const q3 = (trajectory3_->last().degrees_of_freedom().position() -
  //                   reference_position).coordinates().y;
  //Speed const v3 =
  //    trajectory3_->last().degrees_of_freedom().velocity().coordinates().x;
  //std::vector<Vector<Length, EarthMoonOrbitPlane>> const positions3 =
  //    ValuesOf(trajectory3_->Positions(), reference_position);
  //LOG(INFO) << ToMathematicaString(positions3);
  //EXPECT_THAT(positions3.size(), Eq(101));
  //EXPECT_THAT(positions3[25].coordinates().x,
  //            AlmostEquals(0.25 * period_ * v3, 1));
  //EXPECT_THAT(positions3[25].coordinates().y, AlmostEquals(q3, 2));
  //EXPECT_THAT(positions3[50].coordinates().x,
  //            AlmostEquals(0.50 * period_ * v3, 1));
  //EXPECT_THAT(positions3[50].coordinates().y, AlmostEquals(q3, 2));
  //EXPECT_THAT(positions3[75].coordinates().x,
  //            AlmostEquals(0.75 * period_ * v3, 1, 2));
  //EXPECT_THAT(positions3[75].coordinates().y, AlmostEquals(q3, 1));
  //EXPECT_THAT(positions3[100].coordinates().x,
  //            AlmostEquals(1.00 * period_ * v3, 1));
  //EXPECT_THAT(positions3[100].coordinates().y, Eq(q3));
}

TEST_F(EphemerisTest, Sputnik1ToSputnik2) {
  not_null<std::unique_ptr<SolarSystem>> const at_спутник_1_launch =
      SolarSystem::AtСпутник1Launch(
          SolarSystem::Accuracy::kAllBodiesAndOblateness);
  not_null<std::unique_ptr<SolarSystem>> const at_спутник_2_launch =
      SolarSystem::AtСпутник2Launch(
          SolarSystem::Accuracy::kAllBodiesAndOblateness);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies =
      at_спутник_1_launch->massive_bodies();
  std::vector<DegreesOfFreedom<ICRFJ2000Ecliptic>> const initial_state =
      at_спутник_1_launch->initial_state();

  std::vector<not_null<MassiveBody const*>> unowned_bodies;
  for (auto const& body : bodies) {
    unowned_bodies.push_back(body.get());
  }

  Ephemeris<ICRFJ2000Ecliptic>
      ephemeris(
          std::move(bodies),
          initial_state,
          at_спутник_1_launch->time(),
          McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Ecliptic>>(),
          45 * Minute,
          0.1 * Milli(Metre),
          5 * Milli(Metre));

    ephemeris.Prolong(at_спутник_2_launch->time());

  // Upper bounds, tight to the nearest order of magnitude.
  static std::map<SolarSystem::Index, Angle> const expected_angle_error = {{}};
  static std::map<SolarSystem::Index,
                  double> const expected_parent_distance_error = {{}};
  static std::map<SolarSystem::Index,
                  double> const expected_parent_offset_error = {
      {SolarSystem::kAriel, 1E-3},
      {SolarSystem::kDione, 1E-3},
      {SolarSystem::kIo, 1E-3},
      {SolarSystem::kOberon, 1E-3},
      {SolarSystem::kTethys, 1E-3},
      {SolarSystem::kTitania, 1E-3},
      {SolarSystem::kTriton, 1E-4},
      {SolarSystem::kCharon, 1E-4},
      {SolarSystem::kEuropa, 1E-4},
      {SolarSystem::kRhea, 1E-4},
      {SolarSystem::kTitan, 1E-4},
      {SolarSystem::kUmbriel, 1E-4},
      {SolarSystem::kEris, 1E-5},  // NOTE(egg): we may want Dysnomia.
      {SolarSystem::kGanymede, 1E-5},
      {SolarSystem::kIapetus, 1E-5},
      {SolarSystem::kMoon, 1E-5},  // What is this?
      {SolarSystem::kCallisto, 1E-6},
      {SolarSystem::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystem::kPluto, 1E-6},  // NOTE(egg): We are missing Hydra and Nyx.
      {SolarSystem::kVenus, 1E-7},
      {SolarSystem::kEarth, 1E-8},
      {SolarSystem::kJupiter, 1E-8},
      {SolarSystem::kNeptune, 1E-8},
      {SolarSystem::kSaturn, 1E-8},
      {SolarSystem::kUranus, 1E-8},
      {SolarSystem::kMars, 1E-9}};
  static std::map<SolarSystem::Index, double> const expected_position_error = {
      {SolarSystem::kEris, 1E-5},  // NOTE(egg): we may want Dysnomia.
      {SolarSystem::kCharon, 1E-6},
      {SolarSystem::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystem::kPluto, 1E-6},
      {SolarSystem::kTethys, 1E-6},
      {SolarSystem::kAriel, 1E-7},
      {SolarSystem::kDione, 1E-7},
      {SolarSystem::kEuropa, 1E-7},
      {SolarSystem::kIo, 1E-7},
      {SolarSystem::kMoon, 1E-7},
      {SolarSystem::kOberon, 1E-7},
      {SolarSystem::kRhea, 1E-7},
      {SolarSystem::kTitan, 1E-7},
      {SolarSystem::kTitania, 1E-7},
      {SolarSystem::kVenus, 1E-7},
      {SolarSystem::kCallisto, 1E-8},
      {SolarSystem::kEarth, 1E-8},
      {SolarSystem::kGanymede, 1E-8},
      {SolarSystem::kIapetus, 1E-8},
      {SolarSystem::kJupiter, 1E-8},
      {SolarSystem::kNeptune, 1E-8},
      {SolarSystem::kSaturn, 1E-8},
      {SolarSystem::kSun, 1E-8},
      {SolarSystem::kTriton, 1E-8},
      {SolarSystem::kUmbriel, 1E-8},
      {SolarSystem::kUranus, 1E-8},
      {SolarSystem::kMars, 1E-9}};
  static std::map<SolarSystem::Index, double> const expected_velocity_error = {
      {SolarSystem::kAriel, 1E-3},
      {SolarSystem::kCharon, 1E-3},
      {SolarSystem::kDione, 1E-3},
      {SolarSystem::kIo, 1E-3},
      {SolarSystem::kPluto, 1E-3},
      {SolarSystem::kTethys, 1E-3},
      {SolarSystem::kEuropa, 1E-4},
      {SolarSystem::kOberon, 1E-4},
      {SolarSystem::kRhea, 1E-4},
      {SolarSystem::kTitania, 1E-4},
      {SolarSystem::kTriton, 1E-4},
      {SolarSystem::kUmbriel, 1E-4},
      {SolarSystem::kEris, 1E-5},  // NOTE(egg): we may want Dysnomia.
      {SolarSystem::kGanymede, 1E-5},
      {SolarSystem::kTitan, 1E-5},
      {SolarSystem::kUranus, 1E-5},
      {SolarSystem::kCallisto, 1E-6},
      {SolarSystem::kIapetus, 1E-6},
      {SolarSystem::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystem::kMoon, 1E-6},
      {SolarSystem::kSaturn, 1E-6},
      {SolarSystem::kEarth, 1E-7},
      {SolarSystem::kJupiter, 1E-7},
      {SolarSystem::kNeptune, 1E-7},
      {SolarSystem::kSun, 1E-7},
      {SolarSystem::kVenus, 1E-7},
      {SolarSystem::kMars, 1E-8}};

  std::vector<DegreesOfFreedom<ICRFJ2000Ecliptic>> const final_state =
      at_спутник_2_launch->initial_state();
  for (std::size_t i = 0; i < unowned_bodies.size(); ++i) {
    SolarSystem::Index const index = static_cast<SolarSystem::Index>(i);
    ContinuousTrajectory<ICRFJ2000Ecliptic> const& trajectory =
        ephemeris.trajectory(unowned_bodies[i]);
    double const position_error = RelativeError(
        final_state[i].position() - kSolarSystemBarycentre,
        trajectory.EvaluatePosition(at_спутник_2_launch->time(), nullptr) -
            kSolarSystemBarycentre);
    double const velocity_error = RelativeError(
        final_state[i].velocity(),
        trajectory.EvaluateVelocity(at_спутник_2_launch->time(), nullptr));
    EXPECT_THAT(position_error, Lt(expected_position_error.at(index)))
        << SolarSystem::name(i);
    EXPECT_THAT(position_error, Gt(expected_position_error.at(index) / 10.0))
        << SolarSystem::name(i);
    EXPECT_THAT(velocity_error, Lt(expected_velocity_error.at(index)))
        << SolarSystem::name(i);
    EXPECT_THAT(velocity_error, Gt(expected_velocity_error.at(index) / 10.0))
        << SolarSystem::name(i);
    if (i != SolarSystem::kSun) {
      // Look at the error in the position relative to the parent.
      Vector<Length, ICRFJ2000Ecliptic> expected =
          at_спутник_2_launch->trajectories()[i]->
              last().degrees_of_freedom().position() -
          at_спутник_2_launch->trajectories()[SolarSystem::parent(i)]->
              last().degrees_of_freedom().position();
      Vector<Length, ICRFJ2000Ecliptic> actual =
          trajectory.EvaluatePosition(at_спутник_2_launch->time(), nullptr) -
          ephemeris.trajectory(unowned_bodies[SolarSystem::parent(i)]).
              EvaluatePosition(at_спутник_2_launch->time(), nullptr);
      if (expected_angle_error.find(index) != expected_angle_error.end()) {
        Area const product_of_norms = expected.Norm() * actual.Norm();
        Angle const angle = ArcTan(
            Wedge(expected, actual).Norm() / product_of_norms,
            InnerProduct(expected, actual) / product_of_norms);
        EXPECT_THAT(angle / Degree,
                    Gt(expected_angle_error.at(index) / Degree * 0.9))
            << SolarSystem::name(i);
        EXPECT_THAT(angle / Degree,
                    Lt(expected_angle_error.at(index) / Degree * 1.1))
            << SolarSystem::name(i);
      }
      if (expected_parent_distance_error.find(index) !=
          expected_parent_distance_error.end()) {
        double const parent_distance_error = RelativeError(expected.Norm(),
                                                  actual.Norm());
        EXPECT_THAT(parent_distance_error,
                    Lt(expected_parent_distance_error.at(index)))
            << SolarSystem::name(i);
        EXPECT_THAT(parent_distance_error,
                    Gt(expected_parent_distance_error.at(index) / 10.0))
            << SolarSystem::name(i);
      }
      if (expected_parent_offset_error.find(index) !=
          expected_parent_offset_error.end()) {
        double const parent_offset_error =  RelativeError(expected, actual);
        EXPECT_THAT(parent_offset_error,
                    Lt(expected_parent_offset_error.at(index)))
            << SolarSystem::name(i);
        EXPECT_THAT(parent_offset_error,
                    Gt(expected_parent_offset_error.at(index) / 10.0))
            << SolarSystem::name(i);
      }
    }
  }
}

}  // namespace physics
}  // namespace principia
