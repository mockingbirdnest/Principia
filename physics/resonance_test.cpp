#include <map>
#include <vector>

#include "physics/kepler_orbit.hpp"
#include "rigid_motion.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using integrators::McLachlanAtela1992Order5Optimal;
using quantities::astronomy::JulianYear;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using testing_utilities::RelativeError;

namespace physics {

#if !defined(_DEBUG)

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
    elements_[jool_].eccentricity = 0.05;
    elements_[jool_].semimajor_axis = 68'773'560'320 * Metre;
    elements_[jool_].inclination = 1.304 * Degree;
    elements_[jool_].longitude_of_ascending_node = 52 * Degree;
    elements_[jool_].argument_of_periapsis =  0 * Degree;
    elements_[jool_].mean_anomaly = 0.1 * Radian;
    elements_[laythe_].eccentricity = 0;
    elements_[laythe_].semimajor_axis = 27'184'000 * Metre;
    elements_[laythe_].inclination = 0 * Degree;
    elements_[laythe_].longitude_of_ascending_node = 0 * Degree;
    elements_[laythe_].argument_of_periapsis = 0 * Degree;
    elements_[laythe_].mean_anomaly = 3.14 * Radian;
    elements_[vall_].eccentricity = 0;
    elements_[vall_].semimajor_axis = 43'152'000 * Metre;
    elements_[vall_].inclination = 0 * Degree;
    elements_[vall_].longitude_of_ascending_node = 0 * Degree;
    elements_[vall_].argument_of_periapsis = 0 * Degree;
    elements_[vall_].mean_anomaly = 0.9 * Radian;
    elements_[tylo_].eccentricity = 0;
    elements_[tylo_].semimajor_axis = 68'500'000 * Metre;
    elements_[tylo_].inclination = 0.025 * Degree;
    elements_[tylo_].longitude_of_ascending_node = 0 * Degree;
    elements_[tylo_].argument_of_periapsis = 0 * Degree;
    elements_[tylo_].mean_anomaly = 3.14 * Radian;
    elements_[bop_].eccentricity = 0.24;
    elements_[bop_].semimajor_axis = 128'500'000 * Metre;
    elements_[bop_].inclination = 15 * Degree;
    elements_[bop_].longitude_of_ascending_node = 10 * Degree;
    elements_[bop_].argument_of_periapsis = 25 * Degree;
    elements_[bop_].mean_anomaly = 0.9 * Radian;
    elements_[pol_].eccentricity = 0.17;
    elements_[pol_].semimajor_axis = 179'890'000 * Metre;
    elements_[pol_].inclination = 4.25 * Degree;
    elements_[pol_].longitude_of_ascending_node = 2 * Degree;
    elements_[pol_].argument_of_periapsis = 15 * Degree;
    elements_[pol_].mean_anomaly = 0.9 * Radian;
    parents_.emplace(jool_, sun_);
    for (auto const moon : joolian_moons_) {
      parents_.emplace(moon, jool_);
    }
    google::LogToStderr();
  }

  void FixVallMeanAnomaly() {
    // NOTE(egg): In the stock game, this is 0.9 rad.  This is very close to an
    // unstable resonance, and as such no interpretation of the orbital elements
    // will make the system stable for any reasonable length of time.  We make
    // this a stable resonance instead.
    elements_[vall_].mean_anomaly = 0 * Radian;
  }

  // Fills |stock_orbits_|.
  // KSP assumes secondaries are massless.
  void ComputeStockOrbits() {
    for (auto const body : jool_system_) {
      stock_orbits_.emplace(
          body,
          KeplerOrbit<KSP>(*parents_[body],
                           test_particle_,
                           elements_[body],
                           game_epoch_));
    }
  }

  // Keep stock's mean motion when changing the gravitational parameter (instead
  // of keeping the semimajor axis).
  void UseStockMeanMotions() {
    for (auto const body : jool_system_) {
      elements_[body].semimajor_axis = std::experimental::nullopt;
      elements_[body].mean_motion =
          stock_orbits_.at(body).elements_at_epoch().mean_motion;
    }
  }

  DegreesOfFreedom<KSP> StockInitialState(not_null<MassiveBody const*> body) {
    if (body == sun_) {
      return origin_;
    } else {
      return StockInitialState(parents_[body]) +
             stock_orbits_.at(body).StateVectors(game_epoch_);
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
    // Mathematica tends to be slow when dealing with quantities, so we give
    // everything in SI units.
    std::vector<double> times;
    std::vector<std::vector<Vector<double, KSP>>> barycentric_positions;
    for (Instant t = begin; t < end; t += 45 * Minute) {
      auto const position = [&ephemeris, t](
          not_null<MassiveBody const*> body) {
        return ephemeris.trajectory(body)->EvaluatePosition(t, nullptr);
      };

      times.emplace_back((t - game_epoch_) / Second);

      BarycentreCalculator<Position<KSP>, GravitationalParameter>
          jool_system_barycentre;
      for (auto const body : jool_system_) {
        jool_system_barycentre.Add(position(body),
                                   body->gravitational_parameter());
      }
      barycentric_positions.emplace_back();
      for (auto const body : jool_system_) {
        // TODO(egg): when our dynamic frames support that, it would make sense
        // to use a nonrotating dynamic frame centred at the barycentre of the
        // Jool system, instead of computing the barycentre and difference
        // ourselves.
        barycentric_positions.back().emplace_back(
            (position(body) - jool_system_barycentre.Get()) / Metre);
      }
    }
    std::ofstream file;
    file.open(name + "_" + purpose + ".generated.wl");
    file << mathematica::Assign(name + purpose + "q", barycentric_positions);
    file << mathematica::Assign(name + purpose + "t", times);
    file.close();
  }

  // Compute and log the measured periods of the moons.
  void LogPeriods(Ephemeris<KSP> const& ephemeris) {
    auto const position = [this, &ephemeris](
        not_null<MassiveBody const*> body, Instant const& t) {
      return ephemeris.trajectory(body)->EvaluatePosition(t, nullptr);
    };
    auto const barycentre = [this, &position](Instant const& t) {
      return Barycentre<Position<KSP>, Mass>(
          {position(jool_, t), position(laythe_, t), position(vall_, t),
           position(tylo_, t), position(bop_, t), position(pol_, t)},
          {jool_->mass(), laythe_->mass(), vall_->mass(), tylo_->mass(),
           bop_->mass(), pol_->mass()});
    };
    auto const barycentric_position = [this, &barycentre, &position,
                                       &ephemeris](
        not_null<MassiveBody const*> body, Instant const& t) {
      return position(body, t) - barycentre(t);
    };

    LOG(INFO) << barycentric_position(laythe_, game_epoch_);
    LOG(INFO) << barycentric_position(vall_, game_epoch_);
    LOG(INFO) << barycentric_position(tylo_, game_epoch_);
    for (auto const moon : {laythe_, vall_, tylo_}) {
      auto const moon_y = [&barycentric_position, moon,
                           this](Instant const& t) {
        return barycentric_position(moon, t).coordinates().y;
      };

      LOG(INFO) << (moon == laythe_ ? "Laythe" : moon == vall_ ? "Vall"
                                                               : "Tylo");

      Sign const s0(moon_y(game_epoch_));
      Instant t0 = game_epoch_;
      Time const Δt = 45 * Minute;
      while (Sign(moon_y(t0)) == s0) {
        t0 += Δt;
      }
      // The moon crosses the xz plane between t0 and t0 - Δt.
      Instant t1 = t0;
      int const orbits = moon == laythe_ ? 8 : moon == vall_ ? 4 : 2;
      for (int i = 0; i < orbits; ++i) {
        while (Sign(moon_y(t1)) != s0) {
          t1 += Δt;
        }
        // The crossing of the xz plane halfway through the orbit occurs between
        // t1 and t1 - Δt.
        while (Sign(moon_y(t1)) == s0) {
          t1 += Δt;
        }
        // The |i|th orbit ends between t1 and t1 - Δt.
      }
      Time const actual_period =
          (Bisect( moon_y, t1 - Δt, t1) - Bisect(moon_y, t0 - Δt, t0)) / orbits;
      Time const expected_period = (2 * π * Radian) /
                                   *elements_[moon].mean_motion;
      LOG(INFO) << "actual period   : " << actual_period;
      LOG(INFO) << "expected period : " << expected_period;
      LOG(INFO) << "error           :"
                << RelativeError(expected_period, actual_period);
    }
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

  // Interpreting the elements as Jacobi coordinates in the Jool system.
  std::vector<DegreesOfFreedom<KSP>> JacobiInitialStates() {
    // Jool-centric coordinates: a nonrotating inertial frame in which Jool is
    // centred and immobile, for building the Jool system in Jacobi coordinates.
    // We only use this frame at |game_epoch_|.
    using JoolCentric = Frame<serialization::Frame::TestTag,
                              serialization::Frame::TEST1, false>;
    // These coordinate systems have the same axes.
    auto const id = OrthogonalMap<KSP, JoolCentric>::Identity();

        std::map<not_null<MassiveBody const*>, KeplerOrbit<KSP>> orbits;

    orbits.emplace(jool_, KeplerOrbit<KSP>(*sun_,
                                           *jool_,
                                           elements_[jool_],
                                           game_epoch_));



    // The barycentre of the bodies of the Jool system considered so far.
    BarycentreCalculator<DegreesOfFreedom<JoolCentric>, GravitationalParameter>
        inner_system_barycentre;
    // TODO(egg): BarycentreCalculator should just have a method that returns
    // the weight of the whole thing, so we're not accumulating it twice.
    GravitationalParameter inner_system_parameter =
        jool_->gravitational_parameter();
    std::map<not_null<MassiveBody const*>, DegreesOfFreedom<JoolCentric>>
        jool_centric_initial_state;

    // Jool.
    DegreesOfFreedom<JoolCentric> const jool_dof = {JoolCentric::origin,
                                                    Velocity<JoolCentric>()};
    inner_system_barycentre.Add(jool_dof, jool_->gravitational_parameter());
    // The elements of each moon are interpreted as the osculating elements of
    // an orbit around a point mass at the barycentre of Jool and the
    // previously-added moons, so that the state vectors are Jacobi coordinates
    // for the system.
    for (auto const moon : joolian_moons_) {
      jool_centric_initial_state.emplace(
          moon,
          inner_system_barycentre.Get() +
              id(KeplerOrbit<KSP>(MassiveBody(inner_system_parameter),
                                  *moon,
                                  elements_[moon],
                                  game_epoch_).StateVectors(game_epoch_)));
      inner_system_parameter += moon->gravitational_parameter();
      inner_system_barycentre.Add(jool_centric_initial_state.at(moon),
                                  moon->gravitational_parameter());
    }

    // |inner_system_barycentre| is now the barycentre of the whole Jool system.
    // We want that to be placed where dictated by Jool's orbit.
    DegreesOfFreedom<KSP> const jool_barycentre_initial_state =
        origin_ + orbits.at(jool_).StateVectors(game_epoch_);
    // TODO(egg): this is very messy, a constructor from a pair of
    // |DegreesOfFreedom|s like the constructor for |RigidTransformation| from a
    // pair of |Position|s would be nice...
    RigidMotion<JoolCentric, KSP> const to_heliocentric(
        RigidTransformation<JoolCentric, KSP>(
            inner_system_barycentre.Get().position(),
            jool_barycentre_initial_state.position(), id.Inverse()),
        AngularVelocity<JoolCentric>(),
        /*velocity_of_to_frame_origin=*/inner_system_barycentre.Get()
                .velocity() - id(jool_barycentre_initial_state.velocity()));

    // The Sun and Jool.
    std::vector<DegreesOfFreedom<KSP>> initial_states = {
        origin_, to_heliocentric(jool_dof)};
    for (auto const moon : joolian_moons_) {
      initial_states.emplace_back(
          to_heliocentric(jool_centric_initial_state.at(moon)));
    }
    return initial_states;
  }

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> owned_bodies_;
  not_null<MassiveBody const *> const sun_, jool_, laythe_, vall_, tylo_, bop_,
      pol_;
  std::vector<not_null<MassiveBody const*>> bodies_;
  std::vector<not_null<MassiveBody const*>> jool_system_;
  std::vector<not_null<MassiveBody const*>> joolian_moons_;
  std::map<not_null<MassiveBody const*>, KeplerianElements<KSP>> elements_;
  std::map<not_null<MassiveBody const*>, KeplerOrbit<KSP>> stock_orbits_;
  // Nullable second type because we want to use operator[].
  std::map<not_null<MassiveBody const*>, MassiveBody const*> parents_;
  MasslessBody test_particle_;
  // TODO(egg): Frame::unmoving_origin, I have to do this in several places.
  DegreesOfFreedom<KSP> const origin_ = {KSP::origin, Velocity<KSP>()};

  // TODO(egg): this is probably UB, but Point doesn't have constexprs.
  Instant const game_epoch_;
  Instant const reference_ = game_epoch_ + 90 * Day;
  Instant const long_time_ = game_epoch_ + 100 * JulianYear;
  Instant const comparison_ = long_time_ + 90 * Day;

 private:
  not_null<MassiveBody const*> AddBody(GravitationalParameter const& μ) {
    owned_bodies_.emplace_back(make_not_null_unique<MassiveBody>(μ));
    return owned_bodies_.back().get();
  }
};

using ResonanceDeathTest = ResonanceTest;

TEST_F(ResonanceDeathTest, Stock) {
  ComputeStockOrbits();
  UseStockMeanMotions();
  auto ephemeris = MakeEphemeris(StockInitialStates());
  ephemeris.Prolong(reference_);
  LogPeriods(ephemeris);
  LogEphemeris(ephemeris, /*reference=*/true, "stock");
  // Where is thy sting?
  EXPECT_DEATH(
      { ephemeris.Prolong(long_time_); },
      R"regex(Apocalypse occurred at \+8\.22960000000000000e\+06 s)regex");
}

TEST_F(ResonanceTest, Corrected) {
  ComputeStockOrbits();
  UseStockMeanMotions();

  // Instead of putting the moons in a 1:2:4 resonance, put them in a
  // 1:4/φ:16/φ^2 dissonance.
  elements_[vall_].mean_motion =
      *elements_[laythe_].mean_motion / 2.47214;
  *elements_[tylo_].mean_motion =
      *elements_[vall_].mean_motion / 2.47214;

  // Put Bop somewhere further away so it's not kicked out.  A 2:3 mean-motion
  // resonance with Pol works well.
  *elements_[bop_].mean_motion = *elements_[pol_].mean_motion / 1.5;

  auto ephemeris = MakeEphemeris(JacobiInitialStates());
  ephemeris.Prolong(reference_);
  LogPeriods(ephemeris);
  LogEphemeris(ephemeris, /*reference=*/true, "corrected");
  ephemeris.Prolong(long_time_);
  ephemeris.Prolong(comparison_);
  LogEphemeris(ephemeris, /*reference=*/false, "corrected");
}

#endif

}  // namespace physics
}  // namespace principia
