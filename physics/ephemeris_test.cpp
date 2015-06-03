#include "physics/ephemeris.hpp"

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using integrators::McLachlanAtela1992Order5Optimal;
using quantities::Abs;
using quantities::Pow;
using quantities::Sqrt;
using si::Kilogram;
using si::Metre;
using si::Milli;
using si::Second;
using ::testing::Eq;
using ::testing::Lt;

namespace physics {

class EphemerisTest : public testing::Test {
 protected:
  using EarthMoonOrbitPlane = Frame<serialization::Frame::TestTag,
                                    serialization::Frame::TEST, true>;

  void SetUpEarthMoonSystem(
      not_null<std::vector<not_null<std::unique_ptr<MassiveBody>>>*> const
          bodies,
      not_null<std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>>*> const
          initial_state,
      not_null<Position<EarthMoonOrbitPlane>*> const centre_of_mass,
      not_null<Time*> const period) {
    auto earth = std::make_unique<MassiveBody>(6E24 * Kilogram);
    auto moon = std::make_unique<MassiveBody>(7E22 * Kilogram);

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
  std::vector<not_null<std::unique_ptr<MassiveBody>>> bodies;
  std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>> initial_state;
  Position<EarthMoonOrbitPlane> centre_of_mass;
  Time period;
  SetUpEarthMoonSystem(&bodies, &initial_state, &centre_of_mass, &period);

  MassiveBody* const earth = bodies[0].get();
  MassiveBody* const moon = bodies[1].get();

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

}  // namespace physics
}  // namespace principia
