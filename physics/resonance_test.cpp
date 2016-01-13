#include "physics/kepler_orbit.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "physics/solar_system.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using integrators::McLachlanAtela1992Order5Optimal;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;

namespace physics {

class ResonanceTest : public ::testing::Test {};


TEST_F(ResonanceTest, JoolSystem) {
  using KSP =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST, true>;
  MasslessBody test_particle;
  // Gravitational parameters from the KSP wiki.
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  auto const add_body = [&bodies](
      GravitationalParameter const& μ) -> not_null<MassiveBody const*> {
    bodies.emplace_back(make_not_null_unique<MassiveBody>(μ));
    return bodies.back().get();
  };
  auto const sun = add_body(1.1723328E+18 * Pow<3>(Metre) / Pow<2>(Second));
  auto const jool = add_body(2.8252800E+14 * Pow<3>(Metre) / Pow<2>(Second));
  auto const laythe = add_body(1.9620000E+12 * Pow<3>(Metre) / Pow<2>(Second));
  auto const vall = add_body(2.0748150E+11 * Pow<3>(Metre) / Pow<2>(Second));
  auto const tylo = add_body(2.8252800E+12 * Pow<3>(Metre) / Pow<2>(Second));
  auto const bop = add_body(2.4868349E+09 * Pow<3>(Metre) / Pow<2>(Second));
  auto const pol = add_body(7.2170208E+08 * Pow<3>(Metre) / Pow<2>(Second));

  // Elements from the KSP wiki.
  KeplerianElements<KSP> jool_elements;
  jool_elements.eccentricity = 0.05;
  jool_elements.semimajor_axis = 68'773'560'320 * Metre;
  jool_elements.inclination = 1.304 * Degree;
  jool_elements.longitude_of_ascending_node = 52 * Degree;
  jool_elements.argument_of_periapsis =  0 * Degree;
  jool_elements.mean_anomaly = 0.1 * Radian;
  KeplerianElements<KSP> laythe_elements;
  laythe_elements.eccentricity = 0;
  laythe_elements.semimajor_axis = 27'184'000 * Metre;
  laythe_elements.inclination = 0 * Degree;
  laythe_elements.longitude_of_ascending_node = 0 * Degree;
  laythe_elements.argument_of_periapsis = 0 * Degree;
  laythe_elements.mean_anomaly = 3.14 * Radian;
  KeplerianElements<KSP> vall_elements;
  vall_elements.eccentricity = 0;
  vall_elements.semimajor_axis = 43'152'000 * Metre;
  vall_elements.inclination = 0 * Degree;
  vall_elements.longitude_of_ascending_node = 0 * Degree;
  vall_elements.argument_of_periapsis = 0 * Degree;
  vall_elements.mean_anomaly = 0.9 * Radian;
  KeplerianElements<KSP> tylo_elements;
  tylo_elements.eccentricity = 0;
  tylo_elements.semimajor_axis = 68'500'000 * Metre;
  tylo_elements.inclination = 0.025 * Degree;
  tylo_elements.longitude_of_ascending_node = 0 * Degree;
  tylo_elements.argument_of_periapsis = 0 * Degree;
  tylo_elements.mean_anomaly = 3.14 * Radian;
  KeplerianElements<KSP> bop_elements;
  bop_elements.eccentricity = 0.24;
  bop_elements.semimajor_axis = 128'500'000 * Metre;
  bop_elements.inclination = 15 * Degree;
  bop_elements.longitude_of_ascending_node = 10 * Degree;
  bop_elements.argument_of_periapsis = 25 * Degree;
  bop_elements.mean_anomaly = 0.9 * Radian;
  KeplerianElements<KSP> pol_elements;
  pol_elements.eccentricity = 0.17;
  pol_elements.semimajor_axis = 179'890'000 * Metre;
  pol_elements.inclination = 4.25 * Degree;
  pol_elements.longitude_of_ascending_node = 2 * Degree;
  pol_elements.argument_of_periapsis = 15 * Degree;
  pol_elements.mean_anomaly = 0.9 * Radian;

  Instant const game_epoch;
  auto const stock_jool_orbit =
      KeplerOrbit<KSP>(*sun, test_particle, game_epoch, jool_elements);
  auto const stock_laythe_orbit =
      KeplerOrbit<KSP>(*jool, test_particle, game_epoch, laythe_elements);
  auto const stock_vall_orbit =
      KeplerOrbit<KSP>(*jool, test_particle, game_epoch, vall_elements);
  auto const stock_tylo_orbit =
      KeplerOrbit<KSP>(*jool, test_particle, game_epoch, tylo_elements);
  auto const stock_bop_orbit =
      KeplerOrbit<KSP>(*jool, test_particle, game_epoch, bop_elements);
  auto const stock_pol_orbit = 
      KeplerOrbit<KSP>(*jool, test_particle, game_epoch, pol_elements);

  DegreesOfFreedom<KSP> const origin = {KSP::origin, Velocity<KSP>()};
  auto const stock_jool_initial_state =
      origin + stock_jool_orbit.PrimocentricStateVectors(game_epoch);
  Ephemeris<KSP> stock_ephemeris(
      std::move(bodies),
      {origin,
       stock_jool_initial_state,
       stock_jool_initial_state +
           stock_laythe_orbit.PrimocentricStateVectors(game_epoch),
       stock_jool_initial_state +
           stock_vall_orbit.PrimocentricStateVectors(game_epoch),
       stock_jool_initial_state +
           stock_tylo_orbit.PrimocentricStateVectors(game_epoch),
       stock_jool_initial_state +
           stock_bop_orbit.PrimocentricStateVectors(game_epoch),
       stock_jool_initial_state +
           stock_pol_orbit.PrimocentricStateVectors(game_epoch)},
      game_epoch,
      McLachlanAtela1992Order5Optimal<Position<KSP>>(),
      45 * Minute,
      5 * Milli(Metre));
  stock_ephemeris.Prolong(game_epoch + 90 * Day);
  std::vector<Instant> times;
  std::vector<std::vector<Displacement<KSP>>> displacements;
  for (Instant t = game_epoch; t < game_epoch + 90 * Day; t += 45 * Minute) {
    auto const position = [&stock_ephemeris, t](
        not_null<MassiveBody const*> body) {
      return stock_ephemeris.trajectory(body)->EvaluatePosition(t, nullptr);
    };
    auto const barycentre =
        Barycentre<Position<KSP>, Mass>(
            {position(jool), position(laythe), position(vall), position(tylo),
             position(bop), position(pol)},
            {jool->mass(), laythe->mass(), vall->mass(), tylo->mass(),
             bop->mass(), pol->mass()});
    times.emplace_back(t);
    displacements.push_back(
        {position(jool) - barycentre, position(laythe) - barycentre,
         position(vall) - barycentre, position(tylo) - barycentre,
         position(bop) - barycentre, position(pol) - barycentre});
  }
  std::ofstream file;
  file.open("stock_jool.wl");
  file << mathematica::Assign("q", displacements);
  file << mathematica::Assign("t", times);
  file.close();
  // fails.
  stock_ephemeris.Prolong(game_epoch + 100 * Day);
}

}  // namespace physics
}  // namespace principia
