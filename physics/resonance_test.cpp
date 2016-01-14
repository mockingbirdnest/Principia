#include <map>
#include <vector>

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
        pol_(AddBody(7.2170208E+08 * Pow<3>(Metre) / Pow<2>(Second))),
        jool_system_({jool_, laythe_, vall_, tylo_, bop_, pol_}),
        joolian_moons_({laythe_, vall_, tylo_, bop_, pol_}) {
    // Elements from the KSP wiki.
    elements_[jool_].conic.eccentricity = 0.05;
    elements_[jool_].conic.semimajor_axis = 68'773'560'320 * Metre;
    elements_[jool_].inclination = 1.304 * Degree;
    elements_[jool_].longitude_of_ascending_node = 52 * Degree;
    elements_[jool_].argument_of_periapsis =  0 * Degree;
    elements_[jool_].mean_anomaly = 0.1 * Radian;
    elements_[laythe_].conic.eccentricity = 0;
    elements_[laythe_].conic.semimajor_axis = 27'184'000 * Metre;
    elements_[laythe_].inclination = 0 * Degree;
    elements_[laythe_].longitude_of_ascending_node = 0 * Degree;
    elements_[laythe_].argument_of_periapsis = 0 * Degree;
    elements_[laythe_].mean_anomaly = 3.14 * Radian;
    elements_[vall_].conic.eccentricity = 0;
    elements_[vall_].conic.semimajor_axis = 43'152'000 * Metre;
    elements_[vall_].inclination = 0 * Degree;
    elements_[vall_].longitude_of_ascending_node = 0 * Degree;
    elements_[vall_].argument_of_periapsis = 0 * Degree;
    elements_[vall_].mean_anomaly = 0.9 * Radian;
    elements_[tylo_].conic.eccentricity = 0;
    elements_[tylo_].conic.semimajor_axis = 68'500'000 * Metre;
    elements_[tylo_].inclination = 0.025 * Degree;
    elements_[tylo_].longitude_of_ascending_node = 0 * Degree;
    elements_[tylo_].argument_of_periapsis = 0 * Degree;
    elements_[tylo_].mean_anomaly = 3.14 * Radian;
    elements_[bop_].conic.eccentricity = 0.24;
    elements_[bop_].conic.semimajor_axis = 128'500'000 * Metre;
    elements_[bop_].inclination = 15 * Degree;
    elements_[bop_].longitude_of_ascending_node = 10 * Degree;
    elements_[bop_].argument_of_periapsis = 25 * Degree;
    elements_[bop_].mean_anomaly = 0.9 * Radian;
    elements_[pol_].conic.eccentricity = 0.17;
    elements_[pol_].conic.semimajor_axis = 179'890'000 * Metre;
    elements_[pol_].inclination = 4.25 * Degree;
    elements_[pol_].longitude_of_ascending_node = 2 * Degree;
    elements_[pol_].argument_of_periapsis = 15 * Degree;
    elements_[pol_].mean_anomaly = 0.9 * Radian;
    parents_.emplace(jool_, sun_);
    for (auto const moon : joolian_moons_) {
      parents_.emplace(moon, jool_);
    }
    for (auto const body : jool_system_) {
      stock_orbits_.emplace(
          body,
          KeplerOrbit<KSP>(*FindOrDie(parents_, body),
                           test_particle_,
                           game_epoch_,
                           FindOrDie(elements_, body)));
    }
  }

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> owned_bodies_;
  not_null<MassiveBody const *> const sun_, jool_, laythe_, vall_, tylo_, bop_,
      pol_;
  std::vector<not_null<MassiveBody const*>> jool_system_;
  std::vector<not_null<MassiveBody const*>> joolian_moons_;
  std::map<not_null<MassiveBody const*>, KeplerianElements<KSP>> elements_;
  std::map<not_null<MassiveBody const*>, KeplerOrbit<KSP>> stock_orbits_;
  std::map<not_null<MassiveBody const*>, not_null<MassiveBody const*>> parents_;
  MasslessBody test_particle_;
  Instant const game_epoch_;
  DegreesOfFreedom<KSP> const origin_ = {KSP::origin, Velocity<KSP>()};

 private:
  not_null<MassiveBody const*> AddBody(GravitationalParameter const& μ) {
    owned_bodies_.emplace_back(make_not_null_unique<MassiveBody>(μ));
    return owned_bodies_.back().get();
  }
};


TEST_F(ResonanceTest, StockJoolSystem) {
  auto const jool_initial_state =
      origin_ + FindOrDie(stock_orbits_, jool_).PrimocentricStateVectors(game_epoch_);
  Ephemeris<KSP> ephemeris(
      std::move(owned_bodies_),
      {origin_,
       jool_initial_state,
       jool_initial_state +
           FindOrDie(stock_orbits_, laythe_).PrimocentricStateVectors(game_epoch_),
       jool_initial_state +
           FindOrDie(stock_orbits_, vall_).PrimocentricStateVectors(game_epoch_),
       jool_initial_state +
           FindOrDie(stock_orbits_, tylo_).PrimocentricStateVectors(game_epoch_),
       jool_initial_state +
           FindOrDie(stock_orbits_, bop_).PrimocentricStateVectors(game_epoch_),
       jool_initial_state +
           FindOrDie(stock_orbits_, pol_).PrimocentricStateVectors(game_epoch_)},
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
  std::map<not_null<MassiveBody const*>, KeplerOrbit<KSP>> orbits;
  for (auto const body : jool_system_) {
    auto elements = elements_[body];
    elements.conic.semimajor_axis = std::experimental::nullopt;
    elements.conic.mean_motion = FindOrDie(stock_orbits_, body).mean_motion();
    orbits.emplace(body, KeplerOrbit<KSP>(*FindOrDie(parents_, body), *body,
                                          game_epoch_, elements));
  }

  auto const jool_barycentre_initial_state =
      origin_ + FindOrDie(orbits, jool_).PrimocentricStateVectors(game_epoch_);
  auto const laythe_from_barycentre =
      FindOrDie(orbits, laythe_).PrimocentricStateVectors(game_epoch_);
  auto const vall_from_barycentre =
      FindOrDie(orbits, vall_).PrimocentricStateVectors(game_epoch_);
  auto const tylo_from_barycentre =
      FindOrDie(orbits, tylo_).PrimocentricStateVectors(game_epoch_);
  auto const bop_from_barycentre =
      FindOrDie(orbits, bop_).PrimocentricStateVectors(game_epoch_);
  auto const pol_from_barycentre =
      FindOrDie(orbits, pol_).PrimocentricStateVectors(game_epoch_);
  auto const jool_initial_state = jool_barycentre_initial_state -
                                  (laythe_from_barycentre * laythe_->mass() +
                                   vall_from_barycentre * vall_->mass() +
                                   tylo_from_barycentre * tylo_->mass() +
                                   bop_from_barycentre * bop_->mass() +
                                   pol_from_barycentre * pol_->mass()) /
                                      jool_->mass();
  Ephemeris<KSP> ephemeris(
      std::move(owned_bodies_),
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
