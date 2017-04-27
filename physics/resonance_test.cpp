
#include <experimental/filesystem>
#include <map>
#include <string>
#include <vector>

#include "base/file.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "rigid_motion.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using base::Error;
using base::make_not_null_unique;
using base::not_null;
using base::OFStream;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Frame;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Sign;
using geometry::Vector;
using geometry::Velocity;
using integrators::McLachlanAtela1992Order5Optimal;
using numerics::Bisect;
using quantities::GravitationalParameter;
using quantities::Mass;
using quantities::Pow;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::RelativeError;
using ::testing::MatchesRegex;

namespace physics {

class ResonanceTest : public ::testing::Test {
 protected:
  using KSP = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST,
                    /*frame_is_inertial=*/true>;

  ResonanceTest() {
    solar_system_.Initialize(
        SOLUTION_DIR / "astronomy" / "ksp_gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" / "ksp_initial_state_0_0.proto.txt");
    ephemeris_ = solar_system_.MakeEphemeris(
        /*fitting_tolerance=*/5 * Milli(Metre),
        Ephemeris<KSP>::FixedStepParameters(
            McLachlanAtela1992Order5Optimal<Position<KSP>>(),
            /*step=*/45 * Minute));
    kerbol_ = solar_system_.massive_body(*ephemeris_, "Kerbol");
    jool_ = solar_system_.massive_body(*ephemeris_, "Jool");
    pol_ = solar_system_.massive_body(*ephemeris_, "Pol");
    bop_ = solar_system_.massive_body(*ephemeris_, "Bop");
    tylo_ = solar_system_.massive_body(*ephemeris_, "Tylo");
    vall_ = solar_system_.massive_body(*ephemeris_, "Vall");
    laythe_ = solar_system_.massive_body(*ephemeris_, "Laythe");

    bodies_ = {kerbol_, jool_, laythe_, vall_, tylo_, bop_, pol_};
    jool_system_ = {jool_, laythe_, vall_, tylo_, bop_, pol_};
    joolian_moons_ = {laythe_, vall_, tylo_, bop_, pol_};

    parents_.emplace(jool_, kerbol_);
    for (auto const moon : joolian_moons_) {
      parents_.emplace(moon, jool_);
    }

    for (auto const body : bodies_) {
      if (body != kerbol_) {
        elements_[body] = solar_system_.MakeKeplerianElements(
            solar_system_.keplerian_initial_state_message(body->name()).
                elements());
      }
    }

    reference_ = solar_system_.epoch() + 90 * Day;
    long_time_ = solar_system_.epoch() + 100 * JulianYear;
    comparison_ = long_time_ + 90 * Day;

    // This test is mostly a tool for investigating orbit stability, so we want
    // log.
    google::LogToStderr();
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
                           solar_system_.epoch()));
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
    if (body == kerbol_) {
      return origin_;
    } else {
      return StockInitialState(parents_[body]) +
             stock_orbits_.at(body).StateVectors(solar_system_.epoch());
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
                    bool const reference,
                    std::string const& name) {
    Instant const begin = reference ? solar_system_.epoch() : long_time_;
    Instant const end = reference ? reference_ : comparison_;
    std::string const purpose = reference ? "reference" : "comparison";
    // Mathematica tends to be slow when dealing with quantities, so we give
    // everything in SI units.
    std::vector<double> times;
    // Indexed chronologically, then by body.
    std::vector<std::vector<Vector<double, KSP>>> barycentric_positions;
    for (Instant t = begin; t < end; t += 45 * Minute) {
      auto const position = [&ephemeris, t](not_null<MassiveBody const*> body) {
        return ephemeris.trajectory(body)->EvaluatePosition(t);
      };

      times.emplace_back((t - solar_system_.epoch()) / Second);

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
    OFStream file(TEMP_DIR / (name + "_" + purpose + ".generated.wl"));
    file << mathematica::Assign(name + purpose + "q", barycentric_positions);
    file << mathematica::Assign(name + purpose + "t", times);
  }

  // Compute and log the measured periods of the moons.
  void LogPeriods(Ephemeris<KSP> const& ephemeris) {
    auto const position = [this, &ephemeris](
        not_null<MassiveBody const*> body, Instant const& t) {
      return ephemeris.trajectory(body)->EvaluatePosition(t);
    };
    auto const barycentre = [this, &position](Instant const& t) {
      BarycentreCalculator<Position<KSP>, Mass> result;
      for (auto const body : jool_system_) {
        result.Add(position(body, t), body->mass());
      }
      return result.Get();
    };
    auto const barycentric_position =
        [this, &barycentre, &ephemeris, &position](
        not_null<MassiveBody const*> body,
        Instant const& t) {
      return position(body, t) - barycentre(t);
    };

    for (auto const moon : {laythe_, vall_, tylo_}) {
      auto const moon_y =
          [this, &barycentric_position, moon](Instant const& t) {
        return barycentric_position(moon, t).coordinates().y;
      };

      LOG(INFO) << moon->name();

      Sign const s0(moon_y(solar_system_.epoch()));
      Instant t0 = solar_system_.epoch();
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
          (Bisect(moon_y, t1 - Δt, t1) - Bisect(moon_y, t0 - Δt, t0)) / orbits;
      Time const expected_period = (2 * π * Radian) /
                                   *elements_[moon].mean_motion;
      LOG(INFO) << "actual period   : " << actual_period;
      LOG(INFO) << "expected period : " << expected_period;
      LOG(INFO) << "error           : "
                << RelativeError(expected_period, actual_period);
    }
  }

  // Interpreting the elements as Jacobi coordinates in the Jool system.
  std::vector<DegreesOfFreedom<KSP>> JacobiInitialStates() {
    // Jool-centric coordinates: a nonrotating inertial frame in which Jool is
    // centred and immobile, for building the Jool system in Jacobi coordinates.
    // We only use this frame at |solar_system_.epoch()|.
    using JoolCentric = Frame<serialization::Frame::TestTag,
                              serialization::Frame::TEST1,
                              /*is_inertial=*/false>;
    // These coordinate systems have the same axes.
    auto const id = OrthogonalMap<KSP, JoolCentric>::Identity();

        std::map<not_null<MassiveBody const*>, KeplerOrbit<KSP>> orbits;

    orbits.emplace(jool_, KeplerOrbit<KSP>(*kerbol_,
                                           *jool_,
                                           elements_[jool_],
                                           solar_system_.epoch()));

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
                                  solar_system_.epoch())
                     .StateVectors(solar_system_.epoch())));
      inner_system_parameter += moon->gravitational_parameter();
      inner_system_barycentre.Add(jool_centric_initial_state.at(moon),
                                  moon->gravitational_parameter());
    }

    // |inner_system_barycentre| is now the barycentre of the whole Jool system.
    // We want that to be placed where dictated by Jool's orbit.
    DegreesOfFreedom<KSP> const jool_barycentre_initial_state =
        origin_ + orbits.at(jool_).StateVectors(solar_system_.epoch());
    // TODO(egg): this is very messy, a constructor from a pair of
    // |DegreesOfFreedom|s like the constructor for |RigidTransformation| from a
    // pair of |Position|s would be nice...
    RigidMotion<JoolCentric, KSP> const to_heliocentric(
        RigidTransformation<JoolCentric, KSP>(
            inner_system_barycentre.Get().position(),
            jool_barycentre_initial_state.position(),
            id.Inverse()),
        AngularVelocity<JoolCentric>(),
        /*velocity_of_to_frame_origin=*/
            inner_system_barycentre.Get().velocity() -
            id(jool_barycentre_initial_state.velocity()));

    // The Sun and Jool.
    std::vector<DegreesOfFreedom<KSP>> initial_states = {
        origin_, to_heliocentric(jool_dof)};
    for (auto const moon : joolian_moons_) {
      initial_states.emplace_back(
          to_heliocentric(jool_centric_initial_state.at(moon)));
    }
    return initial_states;
  }

  SolarSystem<KSP> solar_system_;
  std::unique_ptr<Ephemeris<KSP>> ephemeris_;

  MassiveBody const* kerbol_;
  MassiveBody const* jool_;
  MassiveBody const* laythe_;
  MassiveBody const* vall_;
  MassiveBody const* tylo_;
  MassiveBody const* bop_;
  MassiveBody const* pol_;
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

  Instant reference_;
  Instant long_time_;
  Instant comparison_;
};

#if !defined(_DEBUG)

TEST_F(ResonanceTest, Stock) {
  ComputeStockOrbits();
  UseStockMeanMotions();
  ephemeris_->Prolong(reference_);
  EXPECT_OK(ephemeris_->last_severe_integration_status());
  LogPeriods(*ephemeris_);
  LogEphemeris(*ephemeris_, /*reference=*/true, "stock");
  ephemeris_->Prolong(long_time_);
  auto const status = ephemeris_->last_severe_integration_status();
  EXPECT_EQ(Error::INVALID_ARGUMENT, status.error());
  EXPECT_THAT(
      status.message(),
      MatchesRegex("Error extending trajectory for Vall\\. Error trying to fit "
                   "a smooth polynomial to the trajectory\\. The approximation "
                   "error jumped from .* m to .* m at time "
                   "\\+8.22960000000000000e\\+06 s\\. The last position is "
                   "\\{.*\\} and the last velocity is \\{.*\\}. An apocalypse "
                   "occurred and two celestials probably collided because your "
                   "solar system is unstable\\."));
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

  ephemeris_->Prolong(reference_);
  LogPeriods(*ephemeris_);
  LogEphemeris(*ephemeris_, /*reference=*/true, "corrected");
  ephemeris_->Prolong(long_time_);
  ephemeris_->Prolong(comparison_);
  LogEphemeris(*ephemeris_, /*reference=*/false, "corrected");
}

#endif

}  // namespace physics
}  // namespace principia
