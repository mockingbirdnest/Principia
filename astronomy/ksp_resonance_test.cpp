
#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "astronomy/stabilize_ksp.hpp"
#include "base/file.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using base::Error;
using base::make_not_null_unique;
using base::not_null;
using base::OFStream;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Frame;
using geometry::Inertial;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Sign;
using geometry::Vector;
using geometry::Velocity;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::McLachlanAtela1992Order5Optimal;
using numerics::Bisect;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::MassiveBody;
using physics::SolarSystem;
using quantities::GravitationalParameter;
using quantities::Infinity;
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
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::operator""_⑴;
using ::testing::Eq;

namespace astronomy {

namespace {
constexpr Time Δt = 45 * Minute;
}  // namespace

class KSPResonanceTest : public ::testing::Test {
 protected:
  using KSP = Frame<enum class KSPTag, KSPTag{}, Inertial>;

  using Periods = std::map<not_null<MassiveBody const*>, Time>;

  KSPResonanceTest()
      : solar_system_(
            SOLUTION_DIR / "astronomy" / "kerbol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" / "kerbol_initial_state_0_0.proto.txt") {
    // This test is mostly a tool for investigating orbit stability, so we want
    // logging.
    google::LogToStderr();
  }

  not_null<std::unique_ptr<Ephemeris<KSP>>> MakeEphemeris() {
    auto ephemeris = solar_system_.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<KSP>::FixedStepParameters(
            SymplecticRungeKuttaNyströmIntegrator<
                McLachlanAtela1992Order5Optimal, Position<KSP>>(),
            /*step=*/Δt));
    jool_ = solar_system_.massive_body(*ephemeris, "Jool");
    laythe_ = solar_system_.massive_body(*ephemeris, "Laythe");
    vall_ = solar_system_.massive_body(*ephemeris, "Vall");
    tylo_ = solar_system_.massive_body(*ephemeris, "Tylo");
    bop_ = solar_system_.massive_body(*ephemeris, "Bop");
    pol_ = solar_system_.massive_body(*ephemeris, "Pol");

    jool_system_ = {jool_, laythe_, vall_, tylo_, bop_, pol_};
    joolian_moons_ = {laythe_, vall_, tylo_, bop_, pol_};

    for (not_null<MassiveBody const*> const moon : joolian_moons_) {
      auto const elements = solar_system_.MakeKeplerianElements(
          solar_system_.keplerian_initial_state_message(moon->name()).
              elements());
      CHECK(elements.mean_motion) << moon->name();
      expected_periods_[moon] = (2 * π * Radian) / *elements.mean_motion;
      longest_joolian_period_ =
          std::max(longest_joolian_period_, expected_periods_[moon]);
    }
    LOG(INFO) << "Longest Joolian period is " << longest_joolian_period_ / Day
              << " days";

    short_term_ = solar_system_.epoch() + 30 * Day;
    mid_term_ = solar_system_.epoch() + 60 * Day;
    long_term_ = solar_system_.epoch() + 100 * JulianYear;

    return ephemeris;
  }

  void LogEphemeris(Ephemeris<KSP> const& ephemeris,
                    Instant const& t_min,
                    Instant const& t_max,
                    std::string const& name) {
    // Mathematica tends to be slow when dealing with quantities, so we give
    // everything in SI units.
    std::vector<double> times;
    // Indexed chronologically, then by body.
    std::vector<std::vector<Vector<double, KSP>>> barycentric_positions;
    for (Instant t = t_min; t <= t_max; t += Δt) {
      auto const position = [&ephemeris, t](not_null<MassiveBody const*> body) {
        return ephemeris.trajectory(body)->EvaluatePosition(t);
      };

      times.emplace_back((t - solar_system_.epoch()) / Second);

      BarycentreCalculator<Position<KSP>, GravitationalParameter>
          jool_system_barycentre;
      for (not_null<MassiveBody const*> const body : jool_system_) {
        jool_system_barycentre.Add(position(body),
                                   body->gravitational_parameter());
      }
      barycentric_positions.emplace_back();
      for (not_null<MassiveBody const*> const body : jool_system_) {
        // TODO(egg): when our dynamic frames support that, it would make sense
        // to use a nonrotating dynamic frame centred at the barycentre of the
        // Jool system, instead of computing the barycentre and difference
        // ourselves.
        barycentric_positions.back().emplace_back(
            (position(body) - jool_system_barycentre.Get()) / Metre);
      }
    }
    OFStream file(TEMP_DIR / (name + ".generated.wl"));
    file << mathematica::Assign(name + "q", barycentric_positions);
    file << mathematica::Assign(name + "t", times);
  }

  // Compute and log the measured periods of the moons.
  Periods ComputePeriods(Ephemeris<KSP> const& ephemeris,
                         Instant const t) const {
    Periods actual_periods;

    auto const position = [&ephemeris](
        not_null<MassiveBody const*> body, Instant const& t) {
      return ephemeris.trajectory(body)->EvaluatePosition(t);
    };
    auto const barycentre = [this, &position](Instant const& t) {
      BarycentreCalculator<Position<KSP>, Mass> result;
      for (not_null<MassiveBody const*> const body : jool_system_) {
        result.Add(position(body, t), body->mass());
      }
      return result.Get();
    };
    auto const barycentric_position =
        [&barycentre, &position](
        not_null<MassiveBody const*> body,
        Instant const& t) {
      return position(body, t) - barycentre(t);
    };

    LOG(INFO) << "Periods at " << t;
    for (auto const moon : {laythe_, vall_, tylo_, pol_, bop_}) {
      auto const moon_y =
          [&barycentric_position, moon](Instant const& t) {
        return barycentric_position(moon, t).coordinates().y;
      };

      Instant t1 = t;
      Sign s0 = Sign::Positive();
      if (t1 <= ephemeris.t_max()) {
        s0 = Sign(moon_y(t1));
      }
      while (t1 <= ephemeris.t_max() && Sign(moon_y(t1)) == s0) {
        t1 += Δt;
      }
      // The moon crosses the xz plane between t1 and t1 - Δt.

      Instant t2 = t1;
      while (t2 <= ephemeris.t_max() && Sign(moon_y(t2)) != s0) {
        t2 += Δt;
      }
      // The crossing of the xz plane halfway through the orbit occurs between
      // t2 and t2 - Δt.
      while (t2 <= ephemeris.t_max() && Sign(moon_y(t2)) == s0) {
        t2 += Δt;
      }
      // The orbit ends between t2 and t2 - Δt.

      LOG(INFO) << "  " << moon->name();
      if (t1 > ephemeris.t_max() || t2 > ephemeris.t_max()) {
        LOG(INFO) << "    Aperiodic";
        actual_periods[moon] = Infinity<Time>();
      } else {
        actual_periods[moon] =
            Bisect(moon_y, t2 - Δt, t2) - Bisect(moon_y, t1 - Δt, t1);
        LOG(INFO) << "    actual period   : " << actual_periods[moon];
        LOG(INFO) << "    expected period : " << expected_periods_.at(moon);
        LOG(INFO) << "    error           : "
                  << RelativeError(actual_periods[moon],
                                   expected_periods_.at(moon));
      }
    }
    return actual_periods;
  }

  SolarSystem<KSP> solar_system_;

  MassiveBody const* jool_;
  MassiveBody const* laythe_;
  MassiveBody const* vall_;
  MassiveBody const* tylo_;
  MassiveBody const* bop_;
  MassiveBody const* pol_;
  std::vector<not_null<MassiveBody const*>> jool_system_;
  std::vector<not_null<MassiveBody const*>> joolian_moons_;
  Time longest_joolian_period_;
  Periods expected_periods_;
  // TODO(egg): Frame::unmoving_origin, I have to do this in several places.
  DegreesOfFreedom<KSP> const origin_ = {KSP::origin, Velocity<KSP>()};

  Instant short_term_;
  Instant mid_term_;
  Instant long_term_;
};

#if !defined(_DEBUG)

TEST_F(KSPResonanceTest, MSVC_ONLY_TEST(Stock)) {
  auto const ephemeris = MakeEphemeris();
  ephemeris->Prolong(short_term_);
  EXPECT_OK(ephemeris->last_severe_integration_status());

  auto const periods_at_epoch =
      ComputePeriods(*ephemeris, ephemeris->t_min());
  EXPECT_THAT(RelativeError(periods_at_epoch.at(laythe_),
                            expected_periods_.at(laythe_)), IsNear(1.4e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_epoch.at(vall_),
                            expected_periods_.at(vall_)), IsNear(2.6e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_epoch.at(tylo_),
                            expected_periods_.at(tylo_)), IsNear(0.91e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_epoch.at(bop_),
                            expected_periods_.at(bop_)), IsNear(9.0e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_epoch.at(pol_),
                            expected_periods_.at(pol_)), IsNear(5.7e-3_⑴));

  auto const periods_at_short_term =
      ComputePeriods(*ephemeris,
                     ephemeris->t_max() - 2 * longest_joolian_period_);
  EXPECT_THAT(RelativeError(periods_at_short_term.at(laythe_),
                            expected_periods_.at(laythe_)), IsNear(1.5e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_short_term.at(vall_),
                            expected_periods_.at(vall_)), IsNear(21e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_short_term.at(tylo_),
                            expected_periods_.at(tylo_)), IsNear(10e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_short_term.at(bop_),
                            expected_periods_.at(bop_)), IsNear(64e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_short_term.at(pol_),
                            expected_periods_.at(pol_)), IsNear(4.2e-3_⑴));

  ephemeris->Prolong(mid_term_);
  EXPECT_OK(ephemeris->last_severe_integration_status());
  auto const periods_at_mid_term =
      ComputePeriods(*ephemeris,
                     ephemeris->t_max() - 2 * longest_joolian_period_);
  EXPECT_THAT(RelativeError(periods_at_mid_term.at(laythe_),
                            expected_periods_.at(laythe_)),
              IsNear(0.874_⑴));
  EXPECT_THAT(periods_at_mid_term.at(vall_), Eq(Infinity<Time>()));
  EXPECT_THAT(periods_at_mid_term.at(tylo_), Eq(Infinity<Time>()));
  EXPECT_THAT(RelativeError(periods_at_mid_term.at(bop_),
                            expected_periods_.at(bop_)), IsNear(0.38_⑴));
  EXPECT_THAT(RelativeError(periods_at_mid_term.at(pol_),
                            expected_periods_.at(pol_)), IsNear(0.30_⑴));

  LogEphemeris(*ephemeris,
               ephemeris->t_max() - 5 * longest_joolian_period_,
               ephemeris->t_max(),
               "stock");
}

TEST_F(KSPResonanceTest, MSVC_ONLY_TEST(Corrected)) {
  StabilizeKSP(solar_system_);

  auto const ephemeris = MakeEphemeris();
  ephemeris->Prolong(short_term_);
  EXPECT_OK(ephemeris->last_severe_integration_status());

  auto const periods_at_epoch =
      ComputePeriods(*ephemeris, ephemeris->t_min());
  EXPECT_THAT(RelativeError(periods_at_epoch.at(laythe_),
                            expected_periods_.at(laythe_)), IsNear(3.9e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_epoch.at(vall_),
                            expected_periods_.at(vall_)), IsNear(5.0e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_epoch.at(tylo_),
                            expected_periods_.at(tylo_)), IsNear(0.74e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_epoch.at(bop_),
                            expected_periods_.at(bop_)), IsNear(13e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_epoch.at(pol_),
                            expected_periods_.at(pol_)), IsNear(11e-3_⑴));

  auto const periods_at_short_term =
      ComputePeriods(*ephemeris,
                     ephemeris->t_max() - 2 * longest_joolian_period_);
  EXPECT_THAT(RelativeError(periods_at_short_term.at(laythe_),
                            expected_periods_.at(laythe_)), IsNear(5.0e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_short_term.at(vall_),
                            expected_periods_.at(vall_)), IsNear(7.8e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_short_term.at(tylo_),
                            expected_periods_.at(tylo_)), IsNear(0.76e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_short_term.at(bop_),
                            expected_periods_.at(bop_)), IsNear(7.6e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_short_term.at(pol_),
                            expected_periods_.at(pol_)), IsNear(22e-3_⑴));

  ephemeris->Prolong(mid_term_);
  EXPECT_OK(ephemeris->last_severe_integration_status());
  auto const periods_at_mid_term =
      ComputePeriods(*ephemeris,
                     ephemeris->t_max() - 2 * longest_joolian_period_);
  EXPECT_THAT(RelativeError(periods_at_mid_term.at(laythe_),
                            expected_periods_.at(laythe_)), IsNear(5.6e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_mid_term.at(vall_),
                            expected_periods_.at(vall_)), IsNear(3.0e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_mid_term.at(tylo_),
                            expected_periods_.at(tylo_)), IsNear(0.73e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_mid_term.at(bop_),
                            expected_periods_.at(bop_)), IsNear(15e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_mid_term.at(pol_),
                            expected_periods_.at(pol_)), IsNear(9.5e-3_⑴));

  ephemeris->Prolong(long_term_);
  EXPECT_OK(ephemeris->last_severe_integration_status());
  auto const periods_at_long_term =
      ComputePeriods(*ephemeris,
                     ephemeris->t_max() - 2 * longest_joolian_period_);
  EXPECT_THAT(RelativeError(periods_at_long_term.at(laythe_),
                            expected_periods_.at(laythe_)), IsNear(4.7e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_long_term.at(vall_),
                            expected_periods_.at(vall_)), IsNear(5.7e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_long_term.at(tylo_),
                            expected_periods_.at(tylo_)), IsNear(0.75e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_long_term.at(bop_),
                            expected_periods_.at(bop_)), IsNear(43e-3_⑴));
  EXPECT_THAT(RelativeError(periods_at_long_term.at(pol_),
                            expected_periods_.at(pol_)), IsNear(9.9e-3_⑴));

  LogEphemeris(*ephemeris,
               ephemeris->t_max() - 5 * longest_joolian_period_,
               ephemeris->t_max(),
               "corrected");
}

#endif

}  // namespace astronomy
}  // namespace principia
