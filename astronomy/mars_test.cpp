#include <fstream>
#include <vector>

#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::Displacement;
using geometry::Position;
using geometry::Velocity;
using physics::ContinuousTrajectory;
using physics::Ephemeris;
using physics::SolarSystem;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;

namespace astronomy {

TEST(MarsTest, Phobos) {
  SolarSystem<ICRFJ2000Equator> solar_system_at_j2000;
  solar_system_at_j2000.Initialize(
      SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "initial_state_jd_2451545_000000000.proto.txt");
  auto const ephemeris =
      solar_system_at_j2000.MakeEphemeris(
          /*fitting_tolerance=*/5 * Milli(Metre),
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              integrators::BlanesMoan2002SRKN14A<Position<ICRFJ2000Equator>>(),
              /*step=*/45 * Minute));
  ephemeris->Prolong(J2000 + 1 * JulianYear);

  ContinuousTrajectory<ICRFJ2000Equator> const& mars_trajectory =
      solar_system_at_j2000.trajectory(*ephemeris, "Mars");
  ContinuousTrajectory<ICRFJ2000Equator>::Hint mars_hint;

  std::vector<Position<ICRFJ2000Equator>> mars_positions;
  std::vector<Velocity<ICRFJ2000Equator>> mars_velocities;

  ContinuousTrajectory<ICRFJ2000Equator> const& phobos_trajectory =
      solar_system_at_j2000.trajectory(*ephemeris, "Phobos");
  ContinuousTrajectory<ICRFJ2000Equator>::Hint phobos_hint;

  std::vector<Position<ICRFJ2000Equator>> phobos_positions;
  std::vector<Velocity<ICRFJ2000Equator>> phobos_velocities;

  std::vector<Displacement<ICRFJ2000Equator>> mars_phobos_displacements;
  std::vector<Velocity<ICRFJ2000Equator>> mars_phobos_velocities;
  for (Instant t = J2000; t < J2000 + 30 * Day; t += 5 * Minute) {
    mars_positions.push_back(mars_trajectory.EvaluatePosition(t, &mars_hint));
    mars_velocities.push_back(mars_trajectory.EvaluateVelocity(t, &mars_hint));
    phobos_positions.push_back(
        phobos_trajectory.EvaluatePosition(t, &phobos_hint));
    phobos_velocities.push_back(
        phobos_trajectory.EvaluateVelocity(t, &phobos_hint));
    mars_phobos_displacements.push_back(
        phobos_positions.back() - mars_positions.back());
    mars_phobos_velocities.push_back(phobos_velocities.back() -
                                     mars_velocities.back());
  }

  std::ofstream file;
  file.open("phobos.generated.wl");
  file << mathematica::Assign("ppaMarsPhobosDisplacements",
                              mars_phobos_displacements);
  file << mathematica::Assign("ppaMarsPhobosVelocities",
                              mars_phobos_velocities);
  file.close();
}

}  // namespace astronomy
}  // namespace principia
