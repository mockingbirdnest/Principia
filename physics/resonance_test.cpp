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

class ResonanceTest : public ::testing::Test {
 protected:
  using KSP =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST, true>;

  // Gravitational parameters from the KSP wiki.
  ResonanceTest()
      : sun_(AddBody(1.1723328E+18 * Pow<3>(Metre) / Pow<2>(Second))),
        jool_(AddBody(2.8252800E+14 * Pow<3>(Metre) / Pow<2>(Second))),
        laythe_(AddBody(1.9620000E+12 * Pow<3>(Metre) / Pow<2>(Second))),
        vall_(AddBody(2.0748150E+11 * Pow<3>(Metre) / Pow<2>(Second))),
        tylo_(AddBody(2.8252800E+12 * Pow<3>(Metre) / Pow<2>(Second))),
        bop_(AddBody(2.4868349E+09 * Pow<3>(Metre) / Pow<2>(Second))),
        pol_(AddBody(7.2170208E+08 * Pow<3>(Metre) / Pow<2>(Second))) {
    // Elements from the KSP wiki.
    jool_elements_.eccentricity = 0.05;
    jool_elements_.semimajor_axis = 68'773'560'320 * Metre;
    jool_elements_.inclination = 1.304 * Degree;
    jool_elements_.longitude_of_ascending_node = 52 * Degree;
    jool_elements_.argument_of_periapsis =  0 * Degree;
    jool_elements_.mean_anomaly = 0.1 * Radian;
    laythe_elements_.eccentricity = 0;
    laythe_elements_.semimajor_axis = 27'184'000 * Metre;
    laythe_elements_.inclination = 0 * Degree;
    laythe_elements_.longitude_of_ascending_node = 0 * Degree;
    laythe_elements_.argument_of_periapsis = 0 * Degree;
    laythe_elements_.mean_anomaly = 3.14 * Radian;
    vall_elements_.eccentricity = 0;
    vall_elements_.semimajor_axis = 43'152'000 * Metre;
    vall_elements_.inclination = 0 * Degree;
    vall_elements_.longitude_of_ascending_node = 0 * Degree;
    vall_elements_.argument_of_periapsis = 0 * Degree;
    vall_elements_.mean_anomaly = 0.9 * Radian;
    tylo_elements_.eccentricity = 0;
    tylo_elements_.semimajor_axis = 68'500'000 * Metre;
    tylo_elements_.inclination = 0.025 * Degree;
    tylo_elements_.longitude_of_ascending_node = 0 * Degree;
    tylo_elements_.argument_of_periapsis = 0 * Degree;
    tylo_elements_.mean_anomaly = 3.14 * Radian;
    bop_elements_.eccentricity = 0.24;
    bop_elements_.semimajor_axis = 128'500'000 * Metre;
    bop_elements_.inclination = 15 * Degree;
    bop_elements_.longitude_of_ascending_node = 10 * Degree;
    bop_elements_.argument_of_periapsis = 25 * Degree;
    bop_elements_.mean_anomaly = 0.9 * Radian;
    pol_elements_.eccentricity = 0.17;
    pol_elements_.semimajor_axis = 179'890'000 * Metre;
    pol_elements_.inclination = 4.25 * Degree;
    pol_elements_.longitude_of_ascending_node = 2 * Degree;
    pol_elements_.argument_of_periapsis = 15 * Degree;
    pol_elements_.mean_anomaly = 0.9 * Radian;
  }

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies_;
  not_null<MassiveBody const *> const sun_, jool_, laythe_, vall_, tylo_, bop_,
      pol_;
  MasslessBody test_particle_;
  Instant const game_epoch_;
  DegreesOfFreedom<KSP> const origin_ = {KSP::origin, Velocity<KSP>()};
  KeplerianElements<KSP> jool_elements_, laythe_elements_, vall_elements_,
      tylo_elements_, bop_elements_, pol_elements_;

 private:
  not_null<MassiveBody const*> AddBody(GravitationalParameter const& μ) {
    bodies_.emplace_back(make_not_null_unique<MassiveBody>(μ));
    return bodies_.back().get();
  }
};


TEST_F(ResonanceTest, JoolSystem) {
  auto const stock_jool_orbit =
      KeplerOrbit<KSP>(*sun_, test_particle_, game_epoch_, jool_elements_);
  auto const stock_laythe_orbit =
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, laythe_elements_);
  auto const stock_vall_orbit =
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, vall_elements_);
  auto const stock_tylo_orbit =
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, tylo_elements_);
  auto const stock_bop_orbit =
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, bop_elements_);
  auto const stock_pol_orbit = 
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, pol_elements_);

  auto const stock_jool_initial_state =
      origin_ + stock_jool_orbit.PrimocentricStateVectors(game_epoch_);
  Ephemeris<KSP> stock_ephemeris(
      std::move(bodies_),
      {origin_,
       stock_jool_initial_state,
       stock_jool_initial_state +
           stock_laythe_orbit.PrimocentricStateVectors(game_epoch_),
       stock_jool_initial_state +
           stock_vall_orbit.PrimocentricStateVectors(game_epoch_),
       stock_jool_initial_state +
           stock_tylo_orbit.PrimocentricStateVectors(game_epoch_),
       stock_jool_initial_state +
           stock_bop_orbit.PrimocentricStateVectors(game_epoch_),
       stock_jool_initial_state +
           stock_pol_orbit.PrimocentricStateVectors(game_epoch_)},
      game_epoch_,
      McLachlanAtela1992Order5Optimal<Position<KSP>>(),
      45 * Minute,
      5 * Milli(Metre));
  stock_ephemeris.Prolong(game_epoch_ + 90 * Day);
  std::vector<Instant> times;
  std::vector<std::vector<Displacement<KSP>>> displacements;
  for (Instant t = game_epoch_; t < game_epoch_ + 90 * Day; t += 45 * Minute) {
    auto const position = [&stock_ephemeris, t](
        not_null<MassiveBody const*> body) {
      return stock_ephemeris.trajectory(body)->EvaluatePosition(t, nullptr);
    };
    auto const barycentre = Barycentre<Position<KSP>, Mass>(
        {position(jool_), position(laythe_), position(vall_), position(tylo_),
         position(bop_), position(pol_)},
        {jool_->mass(), laythe_->mass(), vall_->mass(), tylo_->mass(),
         bop_->mass(), pol_->mass()});
    times.emplace_back(t);
    displacements.push_back(
        {position(jool_) - barycentre, position(laythe_) - barycentre,
         position(vall_) - barycentre, position(tylo_) - barycentre,
         position(bop_) - barycentre, position(pol_) - barycentre});
  }
  std::ofstream file;
  file.open("stock_jool.wl");
  file << mathematica::Assign("q", displacements);
  file << mathematica::Assign("t", times);
  file.close();
  // fails.
  stock_ephemeris.Prolong(game_epoch_ + 100 * Day);
}

}  // namespace physics
}  // namespace principia
