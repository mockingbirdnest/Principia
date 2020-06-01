
#include <algorithm>
#include <chrono>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "absl/strings/ascii.h"
#include "astronomy/stabilize_ksp.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "integrators/integrators.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/hierarchical_system.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/solar_system.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/astronomy.hpp"

namespace principia {

using base::not_null;
using geometry::BarycentreCalculator;
using geometry::Displacement;
using geometry::Frame;
using geometry::Inertial;
using geometry::Instant;
using geometry::Position;
using geometry::Sign;
using geometry::Vector;
using integrators::FixedStepSizeIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order8;
using integrators::methods::QuinlanTremaine1990Order10;
using integrators::methods::QuinlanTremaine1990Order12;
using integrators::methods::BlanesMoan2002SRKN11B;
using integrators::methods::BlanesMoan2002SRKN14A;
using integrators::methods::McLachlanAtela1992Order5Optimal;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::MasslessBody;
using physics::RelativeDegreesOfFreedom;
using physics::SolarSystem;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;
using ::testing::Lt;
using ::testing::Matcher;
using ::testing::_;

namespace astronomy {

using KSP = Frame<enum class KSPTag, Inertial>;

class KSPSystem {
 protected:
  KSPSystem()
      : solar_system_(
            SOLUTION_DIR / "astronomy" / "kerbol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" / "kerbol_initial_state_0_0.proto.txt") {
    StabilizeKSP(solar_system_);
  }

  SolarSystem<KSP> solar_system_;
};

// We apologize for the inheritance.
class KSPSystemTest : public ::testing::Test, protected KSPSystem {
 protected:
  KSPSystemTest()
      : ephemeris_(solar_system_.MakeEphemeris(
            Ephemeris<KSP>::AccuracyParameters(
                /*fitting_tolerance=*/1 * Milli(Metre),
                /*geopotential_tolerance=*/0x1p-24),
            Ephemeris<KSP>::FixedStepParameters(
                SymplecticRungeKuttaNyströmIntegrator<
                    McLachlanAtela1992Order5Optimal,
                    Position<KSP>>(),
                /*step=*/45 * Minute))),
        sun_(solar_system_.massive_body(*ephemeris_, "Sun")),
        eeloo_(solar_system_.massive_body(*ephemeris_, "Eeloo")),
        jool_(solar_system_.massive_body(*ephemeris_, "Jool")),
        pol_(solar_system_.massive_body(*ephemeris_, "Pol")),
        bop_(solar_system_.massive_body(*ephemeris_, "Bop")),
        tylo_(solar_system_.massive_body(*ephemeris_, "Tylo")),
        vall_(solar_system_.massive_body(*ephemeris_, "Vall")),
        laythe_(solar_system_.massive_body(*ephemeris_, "Laythe")),
        dres_(solar_system_.massive_body(*ephemeris_, "Dres")),
        duna_(solar_system_.massive_body(*ephemeris_, "Duna")),
        ike_(solar_system_.massive_body(*ephemeris_, "Ike")),
        kerbin_(solar_system_.massive_body(*ephemeris_, "Kerbin")),
        minmus_(solar_system_.massive_body(*ephemeris_, "Minmus")),
        mun_(solar_system_.massive_body(*ephemeris_, "Mun")),
        eve_(solar_system_.massive_body(*ephemeris_, "Eve")),
        gilly_(solar_system_.massive_body(*ephemeris_, "Gilly")),
        moho_(solar_system_.massive_body(*ephemeris_, "Moho")),
        all_bodies_{sun_,
                    eeloo_,
                    jool_,
                    pol_,
                    bop_,
                    tylo_,
                    vall_,
                    laythe_,
                    dres_,
                    duna_,
                    ike_,
                    kerbin_,
                    minmus_,
                    mun_,
                    eve_,
                    gilly_,
                    moho_},
        planets_and_moons_{eeloo_,
                           jool_,
                           pol_,
                           bop_,
                           tylo_,
                           vall_,
                           laythe_,
                           dres_,
                           duna_,
                           ike_,
                           kerbin_,
                           minmus_,
                           mun_,
                           eve_,
                           gilly_,
                           moho_},
        jool_system_{jool_, laythe_, vall_, tylo_, pol_, bop_},
        joolian_moons_{laythe_, vall_, tylo_, pol_, bop_} {}

  void LogPositions(Ephemeris<KSP> const& ephemeris,
                    Instant const& initial_time,
                    Time const& duration,
                    Matcher<Length> const& matcher,
                    std::string const& name,
                    mathematica::Logger& logger) {
    for (Instant t = initial_time;
         t < initial_time + duration;
         t += 45 * Minute) {
      auto const position =
          [t, &ephemeris](not_null<MassiveBody const*> const body) {
            return ephemeris.trajectory(body)->EvaluatePosition(t);
          };
      BarycentreCalculator<Position<KSP>, GravitationalParameter>
          jool_system_barycentre;
      for (not_null<MassiveBody const*> const body : jool_system_) {
        jool_system_barycentre.Add(position(body),
                                   body->gravitational_parameter());
      }
      std::vector<Displacement<KSP>> barycentric_positions;
      for (not_null<MassiveBody const*> const body : jool_system_) {
        barycentric_positions.emplace_back(
            position(body) - jool_system_barycentre.Get());
        EXPECT_THAT(barycentric_positions.back().Norm(), matcher);
      }
      logger.Append(name, barycentric_positions, mathematica::ExpressIn(Metre));
    }
  }

  not_null<std::unique_ptr<Ephemeris<KSP>>> ephemeris_;

  not_null<MassiveBody const*> const sun_;
  not_null<MassiveBody const*> const eeloo_;
  not_null<MassiveBody const*> const jool_;
  not_null<MassiveBody const*> const pol_;
  not_null<MassiveBody const*> const bop_;
  not_null<MassiveBody const*> const tylo_;
  not_null<MassiveBody const*> const vall_;
  not_null<MassiveBody const*> const laythe_;
  not_null<MassiveBody const*> const dres_;
  not_null<MassiveBody const*> const duna_;
  not_null<MassiveBody const*> const ike_;
  not_null<MassiveBody const*> const kerbin_;
  not_null<MassiveBody const*> const minmus_;
  not_null<MassiveBody const*> const mun_;
  not_null<MassiveBody const*> const eve_;
  not_null<MassiveBody const*> const gilly_;
  not_null<MassiveBody const*> const moho_;
  std::vector<not_null<MassiveBody const*>> const all_bodies_;
  std::vector<not_null<MassiveBody const*>> const planets_and_moons_;
  std::vector<not_null<MassiveBody const*>> const jool_system_;
  std::vector<not_null<MassiveBody const*>> const joolian_moons_;
};

#if !defined(_DEBUG)
TEST_F(KSPSystemTest, KerbalSystem) {
  google::LogToStderr();
  mathematica::Logger logger(TEMP_DIR / "ksp_system.generated.wl",
                             /*make_unique=*/false);

#if 0
  auto const a_century_hence = solar_system_.epoch() + 100 * JulianYear;
#else  // A small century so the tests don't take too long.
  auto const a_century_hence = solar_system_.epoch() + 5 * JulianYear;
#endif

  LOG(INFO) << "Starting integration";
  ephemeris_->Prolong(a_century_hence);
  LOG(INFO) << "Integration done";

  std::map<not_null<MassiveBody const*>, Length> last_separations;
  std::map<not_null<MassiveBody const*>, Sign> last_separation_changes;

  Instant t = solar_system_.epoch();

  for (not_null<MassiveBody const*> const moon : joolian_moons_) {
    last_separation_changes.emplace(moon, Sign::Positive());
  }
  for (int n = 0;
       t < a_century_hence;
       ++n, t = solar_system_.epoch() + n * Hour) {
    auto const position =
        [this, t](not_null<MassiveBody const*> const body) {
          return ephemeris_->trajectory(body)->EvaluatePosition(t);
        };
    auto const degrees_of_freedom =
        [this, t](not_null<MassiveBody const*> const body) {
          return ephemeris_->trajectory(body)->EvaluateDegreesOfFreedom(t);
        };
    auto const jool_position = position(jool_);

    for (not_null<MassiveBody const*> const moon : joolian_moons_) {
      std::string const moon_name = absl::AsciiStrToLower(moon->name());
      Length const separation = (jool_position - position(moon)).Norm();
      Sign const separation_change = Sign(separation - last_separations[moon]);
      if (separation_change != last_separation_changes.at(moon)) {
        logger.Append(moon_name + "Separations",
                      last_separations[moon],
                      mathematica::ExpressIn(Metre));
        logger.Append(moon_name + "Times",
                      t - 1 * Hour - solar_system_.epoch(),
                      mathematica::ExpressIn(Second));
      }
      last_separations[moon] = separation;
      last_separation_changes.at(moon) = separation_change;
    }

    logger.Append("tyloBop",
                  (position(tylo_) - position(bop_)).Norm(),
                  mathematica::ExpressIn(Metre));
    logger.Append("polBop",
                  (position(pol_) - position(bop_)).Norm(),
                  mathematica::ExpressIn(Metre));

    {
      // KSP's osculating elements.
      auto const bop_elements =
          KeplerOrbit<KSP>(*jool_,
                           MasslessBody(),
                           degrees_of_freedom(bop_) - degrees_of_freedom(jool_),
                           t).elements_at_epoch();
      logger.Append("bopEccentricities", *bop_elements.eccentricity);
      logger.Append("bopInclinations",
                    bop_elements.inclination,
                    mathematica::ExpressIn(Degree));
      logger.Append("bopNodes",
                    bop_elements.longitude_of_ascending_node,
                    mathematica::ExpressIn(Degree));
      logger.Append("bopArguments",
                    *bop_elements.argument_of_periapsis,
                    mathematica::ExpressIn(Degree));
    }

    {
      BarycentreCalculator<DegreesOfFreedom<KSP>, GravitationalParameter>
          innermost_jool_system;
      for (not_null<MassiveBody const*> const body :
               {jool_, laythe_, vall_, tylo_}) {
        innermost_jool_system.Add(degrees_of_freedom(body),
                                  body->gravitational_parameter());
      }
      auto const bop_jacobi_elements =
          KeplerOrbit<KSP>(
              MassiveBody(innermost_jool_system.weight()),
              *bop_,
              degrees_of_freedom(bop_) - innermost_jool_system.Get(),
              t).elements_at_epoch();
      logger.Append("bopJacobiEccentricities",
                    *bop_jacobi_elements.eccentricity);
      logger.Append("bopJacobiInclinations",
                    bop_jacobi_elements.inclination,
                    mathematica::ExpressIn(Degree));
      logger.Append("bopJacobiNodes",
                    bop_jacobi_elements.longitude_of_ascending_node,
                    mathematica::ExpressIn(Degree));
      logger.Append("bopJacobiArguments",
                    *bop_jacobi_elements.argument_of_periapsis,
                    mathematica::ExpressIn(Degree));
    }
  }

  LogPositions(*ephemeris_,
               solar_system_.epoch(),
               1 * JulianYear,
               Lt(3e8 * Metre),
               "barycentricPositions1",
               logger);
  LogPositions(*ephemeris_,
                solar_system_.epoch(),
               2 * JulianYear,
               _,
               "barycentricPositions2",
               logger);

}
#endif

struct ConvergenceTestParameters {
  FixedStepSizeIntegrator<
      Ephemeris<KSP>::NewtonianMotionEquation> const& integrator;
  int iterations;
  double first_step_in_seconds;
};

class KSPSystemConvergenceTest
    : public ::testing::TestWithParam<ConvergenceTestParameters>,
      protected KSPSystem {
 public:
  static void SetUpTestCase() {
    logger_ = new mathematica::Logger(
        SOLUTION_DIR / "mathematica" / "ksp_system_convergence.generated.wl",
        /*make_unique=*/false);
  }

  static void TearDownTestCase() {
    delete logger_;
  }

 protected:
  FixedStepSizeIntegrator<
      Ephemeris<KSP>::NewtonianMotionEquation> const&
  integrator() const {
    return GetParam().integrator;
  }

  int iterations() const {
    return GetParam().iterations;
  }

  double first_step_in_seconds() const {
    return GetParam().first_step_in_seconds;
  }

  static mathematica::Logger* logger_;
};

mathematica::Logger* KSPSystemConvergenceTest::logger_ = nullptr;

// This takes 2 minutes to run.
TEST_P(KSPSystemConvergenceTest, DISABLED_Convergence) {
  google::LogToStderr();

  Time const integration_duration = 1 * JulianYear;

  std::map<std::string, std::vector<DegreesOfFreedom<KSP>>>
      name_to_degrees_of_freedom;
  std::vector<Time> steps;
  std::vector<std::chrono::duration<double>> durations;
  for (int i = 0; i < iterations(); ++i) {
    Time const step = first_step_in_seconds() * (1 << i) * Second;
    steps.push_back(step);
    LOG(INFO) << "Integrating with step " << step;

    auto const start = std::chrono::system_clock::now();
    auto const ephemeris = solar_system_.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<KSP>::FixedStepParameters(integrator(), step));
    ephemeris->Prolong(solar_system_.epoch() + integration_duration);
    auto const end = std::chrono::system_clock::now();
    durations.push_back(end - start);

    for (auto const& name : solar_system_.names()) {
      name_to_degrees_of_freedom[name].emplace_back(
          solar_system_.trajectory(*ephemeris, name).
              EvaluateDegreesOfFreedom(solar_system_.epoch() +
                                       integration_duration));
    }
  }

  std::map<std::string, std::vector<RelativeDegreesOfFreedom<KSP>>>
      name_to_errors;
  for (auto const& [name, degrees_of_freedom] : name_to_degrees_of_freedom) {
    CHECK_EQ(degrees_of_freedom.size(), iterations());
    for (int i = 1; i < iterations(); ++i) {
      name_to_errors[name].emplace_back(degrees_of_freedom[i] -
                                        degrees_of_freedom[0]);
    }
  }

  std::vector<Length> position_errors(iterations() - 1);
  std::vector<std::string> worst_body(iterations() - 1);
  std::string const test_name(
      ::testing::UnitTest::GetInstance()->current_test_info()->name());
  for (int i = 0; i < iterations() - 1; ++i) {
    for (auto const& [name, errors] : name_to_errors) {
      if (position_errors[i] < errors[i].displacement().Norm()) {
        worst_body[i] = name;
      }
      position_errors[i] = std::max(position_errors[i],
                                    errors[i].displacement().Norm());
    }
    logger_->Append(std::string("ppaKSPSystemConvergence") +
                        test_name[test_name.size() - 1],
                    std::tuple(steps[i + 1],
                               position_errors[i],
                               worst_body[i],
                               durations[i + 1].count() * Second));
  }
}

INSTANTIATE_TEST_CASE_P(
    AllKSPSystemConvergenceTests,
    KSPSystemConvergenceTest,
    ::testing::Values(
        // This is our preferred integrator.  For a step of 2100 s and an
        // integration over a year, it gives a position error of about 111 m on
        // Laythe and takes about 0.44 s of elapsed time.
        ConvergenceTestParameters{
            SymplecticRungeKuttaNyströmIntegrator<BlanesMoan2002SRKN14A,
                                                  Position<KSP>>(),
            /*iterations=*/7,
            /*first_step_in_seconds=*/65.625},
        ConvergenceTestParameters{
            SymplecticRungeKuttaNyströmIntegrator<
                McLachlanAtela1992Order5Optimal,
                Position<KSP>>(),
            /*iterations=*/8,
            /*first_step_in_seconds=*/32},
        ConvergenceTestParameters{
            SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                               Position<KSP>>(),
            /*iterations=*/6,
            /*first_step_in_seconds=*/64},
        ConvergenceTestParameters{
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order8,
                                               Position<KSP>>(),
            /*iterations=*/6,
            /*first_step_in_seconds=*/64},
        ConvergenceTestParameters{
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order10,
                                               Position<KSP>>(),
            /*iterations=*/6,
            /*first_step_in_seconds=*/64},

        // This is a nice integrator but unfortunately it becomes unstable when
        // Pol and Bop get too close to one another.  For a step of 600 s and an
        // integration over a year, it gives a position error of about 28 m on
        // Bop and takes about 0.7 s of elapsed time.
        ConvergenceTestParameters{
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<KSP>>(),
            /*iterations=*/5,
            /*first_step_in_seconds=*/75}));

}  // namespace astronomy
}  // namespace principia
