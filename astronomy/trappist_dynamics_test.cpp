
#include "astronomy/frames.hpp"
#include "base/file.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::OFStream;
using geometry::Instant;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using physics::Ephemeris;
using physics::KeplerOrbit;
using physics::SolarSystem;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;

namespace astronomy {

class TrappistDynamicsTest : public ::testing::Test {
 protected:
};

TEST_F(TrappistDynamicsTest, MathematicaPeriod) {
  SolarSystem<Trappist> trappist_system(
      SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "trappist_initial_state_jd_2457282_805700000.proto.txt");

  auto const ephemeris = trappist_system.MakeEphemeris(
      /*fitting_tolerance=*/5 * Milli(Metre),
      Ephemeris<Trappist>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                             Position<Trappist>>(),
          /*step=*/0.07 * Day));

  Instant const a_century_later = trappist_system.epoch() + 100 * JulianYear;
  ephemeris->Prolong(a_century_later);

  auto const& star = trappist_system.massive_body(*ephemeris, "Trappist-1A");
  auto const& star_trajectory = ephemeris->trajectory(star);

  OFStream file(TEMP_DIR / "trappist.generated.wl");
  auto const bodies = ephemeris->bodies();
  for (auto const& planet : bodies) {
    if (planet != star) {
      auto const& planet_trajectory = ephemeris->trajectory(planet);
      std::vector<Time> periods;
      for (Instant t = ephemeris->t_max() - 2000 * Hour;
           t < ephemeris->t_max();
           t += 1 * Hour) {
        KeplerOrbit<Trappist> const planet_orbit(
            *star,
            *planet,
            planet_trajectory->EvaluateDegreesOfFreedom(t) -
                star_trajectory->EvaluateDegreesOfFreedom(t),
            t);
        periods.push_back(*planet_orbit.elements_at_epoch().period);
      }

      auto sanitized_name = planet->name();
      sanitized_name.erase(sanitized_name.find_first_of("-"), 1);
      file << mathematica::Assign("period" + sanitized_name,
                                  periods);
    }
  }
}

}  // namespace astronomy
}  // namespace principia
