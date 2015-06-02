#include "physics/ephemeris.hpp"

#include "geometry/frame.hpp"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using integrators::McLachlanAtela1992Order5Optimal;
using quantities::Pow;
using quantities::Sqrt;
using si::Kilogram;
using si::Metre;
using si::Second;

namespace physics {

class EphemerisTest : public testing::Test {
 protected:
  using EarthMoonOrbitPlane = Frame<serialization::Frame::TestTag,
                                    serialization::Frame::TEST, true>;

  EphemerisTest()
      : ephemeris_(std::vector<not_null<std::unique_ptr<MassiveBody>>>(),
                   std::vector<DegreesOfFreedom<EarthMoonOrbitPlane>>(),
                   t0_,
                   McLachlanAtela1992Order5Optimal<Position<EarthMoonOrbitPlane>>(),
                   1 * Second,
                   1 * Metre,
                   10 * Metre) {
  }

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

  Ephemeris<EarthMoonOrbitPlane> ephemeris_;
  Instant t0_;
};

// The canonical Earth-Moon system, tuned to produce circular orbits.
TEST_F(EphemerisTest, EarthMoon) {
  std::vector<Vector<Length, EarthMoonOrbitPlane>> positions;
  system_->Integrate(*integrator_,
                     trajectory1_->last().time() + period_,
                     period_ / 100,
                     1,      // sampling_period
                     false,  // tmax_is_exact
                     {trajectory1_.get(), trajectory2_.get()});

  positions = ValuesOf(trajectory1_->Positions(), centre_of_mass_);
  EXPECT_THAT(positions.size(), Eq(101));
  LOG(INFO) << ToMathematicaString(positions);
  EXPECT_THAT(Abs(positions[25].coordinates().y), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[50].coordinates().x), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[75].coordinates().y), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[100].coordinates().x), Lt(3E-2 * SIUnit<Length>()));

  positions = ValuesOf(trajectory2_->Positions(), centre_of_mass_);
  LOG(INFO) << ToMathematicaString(positions);
  EXPECT_THAT(positions.size(), Eq(101));
  EXPECT_THAT(Abs(positions[25].coordinates().y), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[50].coordinates().x), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[75].coordinates().y), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[100].coordinates().x), Lt(2 * SIUnit<Length>()));
}

TEST_F(EphemerisTest, Test) {
  ephemeris_.Prolong(t0_ + 2 * Second);
}

}  // namespace physics
}  // namespace principia
