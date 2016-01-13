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
    jool_elements_.conic.eccentricity = 0.05;
    jool_elements_.conic.semimajor_axis = 68'773'560'320 * Metre;
    jool_elements_.inclination = 1.304 * Degree;
    jool_elements_.longitude_of_ascending_node = 52 * Degree;
    jool_elements_.argument_of_periapsis =  0 * Degree;
    jool_elements_.mean_anomaly = 0.1 * Radian;
    laythe_elements_.conic.eccentricity = 0;
    laythe_elements_.conic.semimajor_axis = 27'184'000 * Metre;
    laythe_elements_.inclination = 0 * Degree;
    laythe_elements_.longitude_of_ascending_node = 0 * Degree;
    laythe_elements_.argument_of_periapsis = 0 * Degree;
    laythe_elements_.mean_anomaly = 3.14 * Radian;
    vall_elements_.conic.eccentricity = 0;
    vall_elements_.conic.semimajor_axis = 43'152'000 * Metre;
    vall_elements_.inclination = 0 * Degree;
    vall_elements_.longitude_of_ascending_node = 0 * Degree;
    vall_elements_.argument_of_periapsis = 0 * Degree;
    vall_elements_.mean_anomaly = 0.9 * Radian;
    tylo_elements_.conic.eccentricity = 0;
    tylo_elements_.conic.semimajor_axis = 68'500'000 * Metre;
    tylo_elements_.inclination = 0.025 * Degree;
    tylo_elements_.longitude_of_ascending_node = 0 * Degree;
    tylo_elements_.argument_of_periapsis = 0 * Degree;
    tylo_elements_.mean_anomaly = 3.14 * Radian;
    bop_elements_.conic.eccentricity = 0.24;
    bop_elements_.conic.semimajor_axis = 128'500'000 * Metre;
    bop_elements_.inclination = 15 * Degree;
    bop_elements_.longitude_of_ascending_node = 10 * Degree;
    bop_elements_.argument_of_periapsis = 25 * Degree;
    bop_elements_.mean_anomaly = 0.9 * Radian;
    pol_elements_.conic.eccentricity = 0.17;
    pol_elements_.conic.semimajor_axis = 179'890'000 * Metre;
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


TEST_F(ResonanceTest, StockJoolSystem) {
  auto const jool_orbit =
      KeplerOrbit<KSP>(*sun_, test_particle_, game_epoch_, jool_elements_);
  auto const laythe_orbit =
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, laythe_elements_);
  auto const vall_orbit =
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, vall_elements_);
  auto const tylo_orbit =
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, tylo_elements_);
  auto const bop_orbit =
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, bop_elements_);
  auto const pol_orbit = 
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, pol_elements_);

  auto const jool_initial_state =
      origin_ + jool_orbit.PrimocentricStateVectors(game_epoch_);
  Ephemeris<KSP> ephemeris(
      std::move(bodies_),
      {origin_,
       jool_initial_state,
       jool_initial_state + laythe_orbit.PrimocentricStateVectors(game_epoch_),
       jool_initial_state + vall_orbit.PrimocentricStateVectors(game_epoch_),
       jool_initial_state + tylo_orbit.PrimocentricStateVectors(game_epoch_),
       jool_initial_state + bop_orbit.PrimocentricStateVectors(game_epoch_),
       jool_initial_state + pol_orbit.PrimocentricStateVectors(game_epoch_)},
      game_epoch_,
      McLachlanAtela1992Order5Optimal<Position<KSP>>(),
      45 * Minute,
      5 * Milli(Metre));
  ephemeris.Prolong(game_epoch_ + 90 * Day);
  std::vector<Instant> times;
  std::vector<std::vector<Displacement<KSP>>> displacements;
  for (Instant t = game_epoch_; t < game_epoch_ + 90 * Day; t += 45 * Minute) {
    auto const position = [&ephemeris, t](
        not_null<MassiveBody const*> body) {
      return ephemeris.trajectory(body)->EvaluatePosition(t, nullptr);
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
  ephemeris.Prolong(game_epoch_ + 100 * Day);
}

TEST_F(ResonanceTest, BarycentricJoolSystem) {
  auto const jool_orbit =
      KeplerOrbit<KSP>(*sun_, *jool_, game_epoch_, jool_elements_);
  auto const laythe_orbit =
      KeplerOrbit<KSP>(*jool_, test_particle_, game_epoch_, laythe_elements_);
  MassiveBody jool_1(jool_->gravitational_parameter() +
                     laythe_->gravitational_parameter());
  auto const vall_orbit =
      KeplerOrbit<KSP>(jool_1, test_particle_, game_epoch_, vall_elements_);
  MassiveBody jool_2(jool_1.gravitational_parameter() +
                     vall_->gravitational_parameter());
  auto const tylo_orbit =
      KeplerOrbit<KSP>(jool_2, test_particle_, game_epoch_, tylo_elements_);
  MassiveBody jool_3(jool_2.gravitational_parameter() +
                     tylo_->gravitational_parameter());
  auto const bop_orbit =
      KeplerOrbit<KSP>(jool_3, test_particle_, game_epoch_, bop_elements_);
  MassiveBody jool_4(jool_3.gravitational_parameter() +
                     bop_->gravitational_parameter());
  auto const pol_orbit = 
      KeplerOrbit<KSP>(jool_4, test_particle_, game_epoch_, pol_elements_);

  auto const jool_barycentre_initial_state =
      origin_ + jool_orbit.PrimocentricStateVectors(game_epoch_);
  auto const laythe_from_barycentre =
      laythe_orbit.PrimocentricStateVectors(game_epoch_);
  auto const vall_from_barycentre =
      vall_orbit.PrimocentricStateVectors(game_epoch_);
  auto const tylo_from_barycentre =
      tylo_orbit.PrimocentricStateVectors(game_epoch_);
  auto const bop_from_barycentre =
      bop_orbit.PrimocentricStateVectors(game_epoch_);
  auto const pol_from_barycentre =
      pol_orbit.PrimocentricStateVectors(game_epoch_);
  auto const jool_initial_state = jool_barycentre_initial_state -
                                  (laythe_from_barycentre * laythe_->mass() +
                                   vall_from_barycentre * vall_->mass() +
                                   tylo_from_barycentre * tylo_->mass() +
                                   bop_from_barycentre * bop_->mass() +
                                   pol_from_barycentre * pol_->mass()) /
                                      jool_->mass();
  Ephemeris<KSP> ephemeris(
      std::move(bodies_),
      {origin_,  // TODO: not actually here
       jool_initial_state,
       jool_barycentre_initial_state + laythe_from_barycentre,
       jool_barycentre_initial_state + vall_from_barycentre,
       jool_barycentre_initial_state + tylo_from_barycentre,
       jool_barycentre_initial_state + bop_from_barycentre,
       jool_barycentre_initial_state + pol_from_barycentre},
      game_epoch_,
      McLachlanAtela1992Order5Optimal<Position<KSP>>(),
      45 * Minute,
      5 * Milli(Metre));
  ephemeris.Prolong(game_epoch_ + 90 * Day);
  std::vector<Instant> times;
  std::vector<std::vector<Displacement<KSP>>> displacements;
  for (Instant t = game_epoch_; t < game_epoch_ + 90 * Day; t += 45 * Minute) {
    auto const position = [&ephemeris, t](
        not_null<MassiveBody const*> body) {
      return ephemeris.trajectory(body)->EvaluatePosition(t, nullptr);
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
  file.open("corrected_jool.wl");
  file << mathematica::Assign("q", displacements);
  file << mathematica::Assign("t", times);
  file.close();
  // fails.
  ephemeris.Prolong(game_epoch_ + 100 * Day);
}

}  // namespace physics
}  // namespace principia
