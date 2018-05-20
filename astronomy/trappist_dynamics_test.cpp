
#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::Instant;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using physics::Ephemeris;
using physics::SolarSystem;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Metre;
using quantities::si::Milli;

namespace astronomy {

class TrappistDynamicsTest : public ::testing::Test {
 protected:
};

TEST_F(TrappistDynamicsTest, Smoke) {
  SolarSystem<Trappist> trappist_system(
      SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "trappist_initial_state_jd_2457282_805700000.proto.txt");

  auto const ephemeris = trappist_system.MakeEphemeris(
      /*fitting_tolerance=*/5 * Milli(Metre),
      Ephemeris<Trappist>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                             Position<Trappist>>(),
          /*step=*/0.08 * Day));

  Instant const a_century_later = trappist_system.epoch() + 100 * JulianYear;
  ephemeris->Prolong(a_century_later);
}

}  // namespace astronomy
}  // namespace principia
