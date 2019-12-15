
#include <algorithm>
#include <chrono>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "astronomy/stabilize_ksp.hpp"
#include "base/file.hpp"
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
using base::OFStream;
using geometry::BarycentreCalculator;
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

namespace astronomy {

using KSP = Frame<enum class KSPTag, KSPTag{}, Inertial>;

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

  void FillPositions(Ephemeris<KSP> const& ephemeris,
                     Instant const& initial_time,
                     Time const& duration,
                     std::vector<std::vector<Vector<double, KSP>>>& container) {
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
      container.emplace_back();
      for (not_null<MassiveBody const*> const body : jool_system_) {
        container.back().emplace_back(
            (position(body) - jool_system_barycentre.Get()) / Metre);
      }
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

#if 0
  auto const a_century_hence = solar_system_.epoch() + 100 * JulianYear;
#else  // A small century so the tests don't take too long.
  auto const a_century_hence = solar_system_.epoch() + 5 * JulianYear;
#endif

  LOG(INFO) << "Starting integration";
  ephemeris_->Prolong(a_century_hence);
  LOG(INFO) << "Integration done";

  std::map<not_null<MassiveBody const*>, std::vector<double>>
      extremal_separations_in_m;
  std::map<not_null<MassiveBody const*>, std::vector<double>> times_in_s;
  std::map<not_null<MassiveBody const*>, Length> last_separations;
  std::map<not_null<MassiveBody const*>, Sign> last_separation_changes;

  Instant t = solar_system_.epoch();

  // Stock elements.
  std::vector<double> bop_eccentricities;
  std::vector<double> bop_inclinations_in_degrees;
  std::vector<double> bop_nodes_in_degrees;
  std::vector<double> bop_arguments_of_periapsis_in_degrees;

  // Elements around the barycentre of Jool, Laythe, and Vall.
  std::vector<double> bop_jacobi_eccentricities;
  std::vector<double> bop_jacobi_nodes_in_degrees;
  std::vector<double> bop_jacobi_inclinations_in_degrees;
  std::vector<double> bop_jacobi_arguments_of_periapsis_in_degrees;

  std::vector<double> tylo_bop_separations_in_m;
  std::vector<double> pol_bop_separations_in_m;

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
      Length const separation = (jool_position - position(moon)).Norm();
      Sign const separation_change = Sign(separation - last_separations[moon]);
      if (separation_change != last_separation_changes.at(moon)) {
        extremal_separations_in_m[moon].emplace_back(last_separations[moon] /
                                                     Metre);
        times_in_s[moon].emplace_back((t - 1 * Hour - solar_system_.epoch()) /
                                      Second);
      }
      last_separations[moon] = separation;
      last_separation_changes.at(moon) = separation_change;
    }

    tylo_bop_separations_in_m.emplace_back(
        (position(tylo_) - position(bop_)).Norm() / Metre);
    pol_bop_separations_in_m.emplace_back(
        (position(pol_) - position(bop_)).Norm() / Metre);

    {
      // KSP's osculating elements.
      auto const bop_elements =
          KeplerOrbit<KSP>(*jool_,
                           MasslessBody(),
                           degrees_of_freedom(bop_) - degrees_of_freedom(jool_),
                           t).elements_at_epoch();
      bop_eccentricities.emplace_back(*bop_elements.eccentricity);
      bop_inclinations_in_degrees.emplace_back(bop_elements.inclination /
                                               Degree);
      bop_nodes_in_degrees.emplace_back(
          bop_elements.longitude_of_ascending_node / Degree);
      bop_arguments_of_periapsis_in_degrees.emplace_back(
          *bop_elements.argument_of_periapsis / Degree);
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
      bop_jacobi_eccentricities.emplace_back(*bop_jacobi_elements.eccentricity);
      bop_jacobi_inclinations_in_degrees.emplace_back(
          bop_jacobi_elements.inclination / Degree);
      bop_jacobi_nodes_in_degrees.emplace_back(
          bop_jacobi_elements.longitude_of_ascending_node / Degree);
      bop_jacobi_arguments_of_periapsis_in_degrees.emplace_back(
          *bop_jacobi_elements.argument_of_periapsis / Degree);
    }
  }

  std::vector<std::vector<Vector<double, KSP>>> barycentric_positions_1_year;
  FillPositions(*ephemeris_,
                solar_system_.epoch(),
                1 * JulianYear,
                barycentric_positions_1_year);
  std::vector<std::vector<Vector<double, KSP>>> barycentric_positions_2_year;
  FillPositions(*ephemeris_,
                solar_system_.epoch(),
                2 * JulianYear,
                barycentric_positions_2_year);

  for (auto const& body_positions : barycentric_positions_1_year) {
    for (auto const& body_position : body_positions) {
      EXPECT_THAT(body_position.Norm(), Lt(3e8));
    }
  }

  OFStream file(TEMP_DIR / "ksp_system.generated.wl");
  file << mathematica::Assign("laytheTimes", times_in_s[laythe_]);
  file << mathematica::Assign("vallTimes", times_in_s[vall_]);
  file << mathematica::Assign("tyloTimes", times_in_s[tylo_]);
  file << mathematica::Assign("polTimes", times_in_s[pol_]);
  file << mathematica::Assign("bopTimes", times_in_s[bop_]);
  file << mathematica::Assign("laytheSeparations",
                              extremal_separations_in_m[laythe_]);
  file << mathematica::Assign("vallSeparations",
                              extremal_separations_in_m[vall_]);
  file << mathematica::Assign("tyloSeparations",
                              extremal_separations_in_m[tylo_]);
  file << mathematica::Assign("polSeparations",
                              extremal_separations_in_m[pol_]);
  file << mathematica::Assign("bopSeparations",
                              extremal_separations_in_m[bop_]);

  file << mathematica::Assign("bopEccentricities", bop_eccentricities);
  file << mathematica::Assign("bopInclinations", bop_inclinations_in_degrees);
  file << mathematica::Assign("bopNodes", bop_nodes_in_degrees);
  file << mathematica::Assign("bopArguments",
                              bop_arguments_of_periapsis_in_degrees);
  file << mathematica::Assign("bopJacobiEccentricities",
                              bop_jacobi_eccentricities);
  file << mathematica::Assign("bopJacobiInclinations",
                              bop_jacobi_inclinations_in_degrees);
  file << mathematica::Assign("bopJacobiNodes", bop_jacobi_nodes_in_degrees);
  file << mathematica::Assign("bopJacobiArguments",
                              bop_jacobi_arguments_of_periapsis_in_degrees);

  file << mathematica::Assign("tyloBop", tylo_bop_separations_in_m);
  file << mathematica::Assign("polBop", pol_bop_separations_in_m);

  file << mathematica::Assign("barycentricPositions1",
                              barycentric_positions_1_year);
  file << mathematica::Assign("barycentricPositions2",
                              barycentric_positions_2_year);
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
    file_ = OFStream(SOLUTION_DIR / "mathematica" /
                     "ksp_system_convergence.generated.wl");
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

  static OFStream file_;
};

OFStream KSPSystemConvergenceTest::file_;

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

  using MathematicaEntry = std::tuple<Time, Length, std::string, Time>;
  using MathematicaEntries = std::vector<MathematicaEntry>;

  std::vector<Length> position_errors(iterations() - 1);
  std::vector<std::string> worst_body(iterations() - 1);
  MathematicaEntries mathematica_entries;
  for (int i = 0; i < iterations() - 1; ++i) {
    for (auto const& [name, errors] : name_to_errors) {
      if (position_errors[i] < errors[i].displacement().Norm()) {
        worst_body[i] = name;
      }
      position_errors[i] = std::max(position_errors[i],
                                    errors[i].displacement().Norm());
    }
    mathematica_entries.push_back({steps[i + 1],
                                   position_errors[i],
                                   mathematica::Escape(worst_body[i]),
                                   durations[i + 1].count() * Second});
  }

  std::string const test_name(
      ::testing::UnitTest::GetInstance()->current_test_info()->name());
  file_ << mathematica::Assign(
      std::string("ppaKSPSystemConvergence") + test_name[test_name.size() - 1],
      mathematica::ToMathematica(mathematica_entries));
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
