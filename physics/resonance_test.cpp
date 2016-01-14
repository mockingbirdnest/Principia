#include <map>
#include <vector>

#include "physics/kepler_orbit.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using integrators::McLachlanAtela1992Order5Optimal;
using quantities::astronomy::JulianYear;
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
        bodies_({sun_, jool_, laythe_, vall_, tylo_, bop_, pol_}),
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
  }

  void FixVallMeanAnomaly() {
    // NOTE(egg): In the stock game, this is 0.9 rad.  This is very close to an
    // unstable resonance, and as such no interpretation of the orbital elements
    // will make the system stable for any reasonable length of time.  We make
    // this a stable resonance instead.
    elements_[vall_].mean_anomaly = 0 * Radian;
  }

  // KSP assumes secondaries are massless.
  void ComputeStockOrbits() {
    for (auto const body : jool_system_) {
      stock_orbits_.emplace(
          body,
          KeplerOrbit<KSP>(*FindOrDie(parents_, body),
                           test_particle_,
                           game_epoch_,
                           FindOrDie(elements_, body)));
    }
  }

  DegreesOfFreedom<KSP> StockInitialState(not_null<MassiveBody const*> body) {
    if (body == sun_) {
      return origin_;
    } else {
      return StockInitialState(parents_.at(body)) +
             stock_orbits_.at(body).PrimocentricStateVectors(game_epoch_);
    }
  }

  std::vector<DegreesOfFreedom<KSP>> StockInitialStates() {
    std::vector<DegreesOfFreedom<KSP>> initial_states;
    for (auto const body : bodies_) {
      initial_states.emplace_back(StockInitialState(body));
    }
    return initial_states;
  }

  void LogEphemeris(Ephemeris<KSP> const& ephemeris,
                    bool reference,
                    std::string name) {
    Instant const begin = reference ? game_epoch_ : long_time_;
    Instant const end = reference ? reference_ : comparison_;
    std::string const purpose = reference ? "reference" : "comparison";
    std::vector<Instant> times;
    std::vector<std::vector<Displacement<KSP>>> displacements;
    std::vector<std::vector<Vector<double, KSP>>> unitless_displacements;
    for (Instant t = game_epoch_; t < reference_; t += 45 * Minute) {
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
      unitless_displacements.emplace_back();
      unitless_displacements.back().resize(displacements.back().size());
      std::transform(displacements.back().begin(), displacements.back().end(),
                     unitless_displacements.back().begin(),
                     [](Displacement<KSP> d) { return d / Metre; });
    }
    std::ofstream file;
    file.open(name + "_" + purpose + ".wl");
    file << mathematica::Assign("q", displacements);
    file << mathematica::Assign("qSI", unitless_displacements);
    file << mathematica::Assign("t", times);
    file.close();
  }

  Ephemeris<KSP> MakeEphemeris(std::vector<DegreesOfFreedom<KSP>> states) {
    return Ephemeris<KSP>(
      std::move(owned_bodies_),
      states,
      game_epoch_,
      McLachlanAtela1992Order5Optimal<Position<KSP>>(),
      45 * Minute,
      5 * Milli(Metre));
  }

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> owned_bodies_;
  not_null<MassiveBody const *> const sun_, jool_, laythe_, vall_, tylo_, bop_,
      pol_;
  std::vector<not_null<MassiveBody const*>> bodies_;
  std::vector<not_null<MassiveBody const*>> jool_system_;
  std::vector<not_null<MassiveBody const*>> joolian_moons_;
  std::map<not_null<MassiveBody const*>, KeplerianElements<KSP>> elements_;
  std::map<not_null<MassiveBody const*>, KeplerOrbit<KSP>> stock_orbits_;
  std::map<not_null<MassiveBody const*>, not_null<MassiveBody const*>> parents_;
  MasslessBody test_particle_;
  DegreesOfFreedom<KSP> const origin_ = {KSP::origin, Velocity<KSP>()};

  // TODO(egg): this is probably UB, but Point doesn't have constexprs.
  Instant const game_epoch_;
  Instant const reference_ = game_epoch_ + 30 * Day;
  Instant const long_time_ = game_epoch_ + 100 * JulianYear;
  Instant const comparison_ = long_time_ + 30 * Day;

 private:
  not_null<MassiveBody const*> AddBody(GravitationalParameter const& μ) {
    owned_bodies_.emplace_back(make_not_null_unique<MassiveBody>(μ));
    return owned_bodies_.back().get();
  }
};


TEST_F(ResonanceTest, StockJoolSystem) {
  ComputeStockOrbits();
  auto ephemeris = MakeEphemeris(StockInitialStates());
  ephemeris.Prolong(reference_);
  LogEphemeris(ephemeris, /*reference=*/true, "stock_jool");
  // where is thy sting
  ephemeris.Prolong(long_time_);
}

TEST_F(ResonanceTest, FixedVallAnomalyJoolSystem) {
  FixVallMeanAnomaly();
  ComputeStockOrbits();
  auto ephemeris = MakeEphemeris(StockInitialStates());
  ephemeris.Prolong(reference_);
  LogEphemeris(ephemeris, /*reference=*/true, "fixed_vall_anomaly_jool");
  ephemeris.Prolong(long_time_);
  ephemeris.Prolong(comparison_);
  LogEphemeris(ephemeris, /*reference=*/false, "fixed_vall_anomaly_jool");
}

TEST_F(ResonanceTest, BarycentricJoolSystem) {
  FixVallMeanAnomaly();
  ComputeStockOrbits();

  std::map<not_null<MassiveBody const*>, RelativeDegreesOfFreedom<KSP>>
      jooliocentric_initial_states;

  KeplerOrbit<KSP> jool_orbit(*sun_, *jool_, game_epoch_, elements_[jool_]);

  GravitationalParameter inner_system_parameter;
  inner_system_parameter = jool_->gravitational_parameter();
  BarycentreCalculator<RelativeDegreesOfFreedom<KSP>, GravitationalParameter>
      inner_system_jooliocentric_barycentre;
  inner_system_jooliocentric_barycentre.Add(
      origin_ - origin_,  // why is this not default-constructible?
      jool_->gravitational_parameter());

  for (auto const body : joolian_moons_) {
    auto elements = elements_[body];
    elements.conic.semimajor_axis = std::experimental::nullopt;
    elements.conic.mean_motion = FindOrDie(stock_orbits_, body).mean_motion();
    LOG(ERROR) << *elements.conic.mean_motion;
    jooliocentric_initial_states.emplace(
        body, inner_system_jooliocentric_barycentre.Get() +
                  KeplerOrbit<KSP>(MassiveBody(inner_system_parameter), *body,
                                   game_epoch_, elements)
                      .PrimocentricStateVectors(game_epoch_));
    inner_system_parameter += body->gravitational_parameter();
    inner_system_jooliocentric_barycentre.Add(
        jooliocentric_initial_states.at(body), body->gravitational_parameter());
  }

  auto const jool_system_jooliocentric_barycentre =
      inner_system_jooliocentric_barycentre.Get();

  auto const jool_barycentre_initial_state =
      origin_ + jool_orbit.PrimocentricStateVectors(game_epoch_);
  auto const jool_initial_state =
      jool_barycentre_initial_state - jool_system_jooliocentric_barycentre;

  auto ephemeris =
      MakeEphemeris(
          {origin_,
           jool_initial_state,
           jool_initial_state + jooliocentric_initial_states.at(laythe_),
           jool_initial_state + jooliocentric_initial_states.at(vall_),
           jool_initial_state + jooliocentric_initial_states.at(tylo_),
           jool_initial_state + jooliocentric_initial_states.at(bop_),
           jool_initial_state + jooliocentric_initial_states.at(pol_)});

  ephemeris.Prolong(reference_);
  LogEphemeris(ephemeris, /*reference=*/true, "barycentric_jool");
  ephemeris.Prolong(long_time_);
  ephemeris.Prolong(comparison_);
  LogEphemeris(ephemeris, /*reference=*/false, "barycentric_jool");
}

}  // namespace physics
}  // namespace principia
