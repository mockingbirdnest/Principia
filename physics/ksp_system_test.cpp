
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/hierarchical_system.hpp"
#include "physics/kepler_orbit.hpp"
#include "quantities/astronomy.hpp"
#include "rigid_motion.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using quantities::astronomy::JulianYear;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using testing_utilities::Componentwise;
using testing_utilities::RelativeError;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;

namespace physics {

class KSPSystemTest : public ::testing::Test {
 protected:
  using KSP = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST,
                    /*frame_is_inertial=*/true>;

  struct KSPCelestial {
    KeplerianElements<KSP> elements;
    KSPCelestial* parent = nullptr;
    MassiveBody* body;
    std::unique_ptr<MassiveBody> owned_body;
  };

  KSPSystemTest() {
    sun_.owned_body = std::make_unique<MassiveBody>(
        1.1723327948324908E+18 * SIUnit<GravitationalParameter>());
    eeloo_.owned_body = std::make_unique<MassiveBody>(
        74410814527.049576 * SIUnit<GravitationalParameter>());
    jool_.owned_body = std::make_unique<MassiveBody>(
        282528004209995.31 * SIUnit<GravitationalParameter>());
    pol_.owned_body = std::make_unique<MassiveBody>(
        721702080.00000012 * SIUnit<GravitationalParameter>());
    bop_.owned_body = std::make_unique<MassiveBody>(
        2486834944.414907 * SIUnit<GravitationalParameter>());
    tylo_.owned_body = std::make_unique<MassiveBody>(
        2825280042099.9531 * SIUnit<GravitationalParameter>());
    vall_.owned_body = std::make_unique<MassiveBody>(
        207481499473.75098 * SIUnit<GravitationalParameter>());
    laythe_.owned_body = std::make_unique<MassiveBody>(
        1962000029236.0784 * SIUnit<GravitationalParameter>());
    dres_.owned_body = std::make_unique<MassiveBody>(
        21484488600.000004 * SIUnit<GravitationalParameter>());
    duna_.owned_body = std::make_unique<MassiveBody>(
        301363211975.09772 * SIUnit<GravitationalParameter>());
    ike_.owned_body = std::make_unique<MassiveBody>(
        18568368573.144012 * SIUnit<GravitationalParameter>());
    kerbin_.owned_body = std::make_unique<MassiveBody>(
        3531600000000 * SIUnit<GravitationalParameter>());
    minmus_.owned_body = std::make_unique<MassiveBody>(
        1765800026.3124719 * SIUnit<GravitationalParameter>());
    mun_.owned_body = std::make_unique<MassiveBody>(
        65138397520.780701 * SIUnit<GravitationalParameter>());
    eve_.owned_body = std::make_unique<MassiveBody>(
        8171730229210.874 * SIUnit<GravitationalParameter>());
    gilly_.owned_body = std::make_unique<MassiveBody>(
        8289449.814716354 * SIUnit<GravitationalParameter>());
    moho_.owned_body = std::make_unique<MassiveBody>(
        168609378654.50949 * SIUnit<GravitationalParameter>());
    for (auto const celestial : all_bodies_) {
      celestial->body = celestial->owned_body.get();
    }
    eeloo_.parent = &sun_;
    jool_.parent = &sun_;
    pol_.parent = &jool_;
    bop_.parent = &jool_;
    tylo_.parent = &jool_;
    vall_.parent = &jool_;
    laythe_.parent = &jool_;
    dres_.parent = &sun_;
    duna_.parent = &sun_;
    ike_.parent = &duna_;
    kerbin_.parent = &sun_;
    minmus_.parent = &kerbin_;
    mun_.parent = &kerbin_;
    eve_.parent = &sun_;
    gilly_.parent = &eve_;
    moho_.parent = &sun_;
    eeloo_.elements.eccentricity = +2.60000000000000009e-01;
    eeloo_.elements.mean_motion = +4.00223155970064009e-08 * (Radian / Second);
    eeloo_.elements.inclination = +1.07337748997651278e-01 * Radian;
    eeloo_.elements.longitude_of_ascending_node =
        +8.72664625997164767e-01 * Radian;
    eeloo_.elements.argument_of_periapsis = +4.53785605518525692e+00 * Radian;
    eeloo_.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    jool_.elements.eccentricity = +5.00000007450581013e-02;
    jool_.elements.mean_motion = +6.00334352457231946e-08 * (Radian / Second);
    jool_.elements.inclination = +2.27590937955459392e-02 * Radian;
    jool_.elements.longitude_of_ascending_node =
        +9.07571211037051406e-01 * Radian;
    jool_.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    jool_.elements.mean_anomaly = +1.00000001490115994e-01 * Radian;
    pol_.elements.eccentricity = +1.70850000000000002e-01;
    pol_.elements.mean_motion = +6.96658945572122982e-06 * (Radian / Second);
    pol_.elements.inclination = +7.41764932097590118e-02 * Radian;
    pol_.elements.longitude_of_ascending_node =
        +3.49065850398865909e-02 * Radian;
    pol_.elements.argument_of_periapsis = +2.61799387799149408e-01 * Radian;
    pol_.elements.mean_anomaly = +8.99999976158141979e-01 * Radian;
    bop_.elements.eccentricity = +2.34999999403953996e-01;
    bop_.elements.mean_motion = +9.95227065103033049e-06 * (Radian / Second);
    bop_.elements.inclination = +2.87979326579064354e+00 * Radian;
    bop_.elements.longitude_of_ascending_node =
        +1.74532925199432948e-01 * Radian;
    bop_.elements.argument_of_periapsis = +4.36332312998582383e-01 * Radian;
    bop_.elements.mean_anomaly = +8.99999976158141979e-01 * Radian;
    tylo_.elements.eccentricity = +0.00000000000000000e+00;
    tylo_.elements.mean_motion = +1.94051054171045988e-05 * (Radian / Second);
    tylo_.elements.inclination = +4.36332319500439990e-04 * Radian;
    tylo_.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    tylo_.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    tylo_.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    vall_.elements.eccentricity = +0.00000000000000000e+00;
    vall_.elements.mean_motion = +4.79720588121814983e-05 * (Radian / Second);
    vall_.elements.inclination = +0.00000000000000000e+00 * Radian;
    vall_.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    vall_.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    vall_.elements.mean_anomaly = +8.99999976158141979e-01 * Radian;
    laythe_.elements.eccentricity = +0.00000000000000000e+00;
    laythe_.elements.mean_motion = +1.18593451424947995e-04 * (Radian / Second);
    laythe_.elements.inclination = +0.00000000000000000e+00 * Radian;
    laythe_.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    laythe_.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    laythe_.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    dres_.elements.eccentricity = +1.44999999999999990e-01;
    dres_.elements.mean_motion = +1.31191970097993002e-07 * (Radian / Second);
    dres_.elements.inclination = +8.72664625997164739e-02 * Radian;
    dres_.elements.longitude_of_ascending_node =
        +4.88692190558412243e+00 * Radian;
    dres_.elements.argument_of_periapsis = +1.57079632679489656e+00 * Radian;
    dres_.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    duna_.elements.eccentricity = +5.09999990463256975e-02;
    duna_.elements.mean_motion = +3.62866884706430976e-07 * (Radian / Second);
    duna_.elements.inclination = +1.04719752778990849e-03 * Radian;
    duna_.elements.longitude_of_ascending_node =
        +2.36492113645231639e+00 * Radian;
    duna_.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    duna_.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    ike_.elements.eccentricity = +2.99999993294477012e-02;
    ike_.elements.mean_motion = +9.59003407994517016e-05 * (Radian / Second);
    ike_.elements.inclination = +3.49065855600351992e-03 * Radian;
    ike_.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    ike_.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    ike_.elements.mean_anomaly = +1.70000004768372004e+00 * Radian;
    kerbin_.elements.eccentricity = +0.00000000000000000e+00;
    kerbin_.elements.mean_motion = +6.82691894080843017e-07 * (Radian / Second);
    kerbin_.elements.inclination = +0.00000000000000000e+00 * Radian;
    kerbin_.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    kerbin_.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    kerbin_.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    minmus_.elements.eccentricity = +0.00000000000000000e+00;
    minmus_.elements.mean_motion = +5.83228807719951003e-06 * (Radian / Second);
    minmus_.elements.inclination = +1.04719755119659780e-01 * Radian;
    minmus_.elements.longitude_of_ascending_node =
        +1.36135681655557694e+00 * Radian;
    minmus_.elements.argument_of_periapsis = +6.63225115757845263e-01 * Radian;
    minmus_.elements.mean_anomaly = +8.99999976158141979e-01 * Radian;
    mun_.elements.eccentricity = +0.00000000000000000e+00;
    mun_.elements.mean_motion = +4.52078533000627999e-05 * (Radian / Second);
    mun_.elements.inclination = +0.00000000000000000e+00 * Radian;
    mun_.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    mun_.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    mun_.elements.mean_anomaly = +1.70000004768372004e+00 * Radian;
    eve_.elements.eccentricity = +9.99999977648258036e-03;
    eve_.elements.mean_motion = +1.11049676511037010e-06 * (Radian / Second);
    eve_.elements.inclination = +3.66519126274052684e-02 * Radian;
    eve_.elements.longitude_of_ascending_node =
        +2.61799387799149408e-01 * Radian;
    eve_.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    eve_.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    gilly_.elements.eccentricity = +5.50000011920928955e-01;
    gilly_.elements.mean_motion = +1.61692985452753988e-05 * (Radian / Second);
    gilly_.elements.inclination = +2.09439510239319560e-01 * Radian;
    gilly_.elements.longitude_of_ascending_node =
        +1.39626340159546358e+00 * Radian;
    gilly_.elements.argument_of_periapsis = +1.74532925199432948e-01 * Radian;
    gilly_.elements.mean_anomaly = +8.99999976158141979e-01 * Radian;
    moho_.elements.eccentricity = +2.00000002980231989e-01;
    moho_.elements.mean_motion = +2.83568694188237007e-06 * (Radian / Second);
    moho_.elements.inclination = +1.22173047639603072e-01 * Radian;
    moho_.elements.longitude_of_ascending_node =
        +1.22173047639603061e+00 * Radian;
    moho_.elements.argument_of_periapsis = +2.61799387799149408e-01 * Radian;
    moho_.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
  }

  not_null<std::unique_ptr<Ephemeris<KSP>>> MakeEphemeris() {
    HierarchicalSystem<KSP> hierarchical_system(std::move(sun_.owned_body));
    for (auto const celestial : planets_and_moons_) {
      hierarchical_system.Add(std::move(celestial->owned_body),
                              celestial->parent->body,
                              celestial->elements);
    }
    HierarchicalSystem<KSP>::BarycentricSystem barycentric_system =
        hierarchical_system.ConsumeBarycentricSystem();
    return make_not_null_unique<Ephemeris<KSP>>(
        std::move(barycentric_system.bodies),
        std::move(barycentric_system.degrees_of_freedom),
        ksp_epoch,
        integrators::McLachlanAtela1992Order5Optimal<Position<KSP>>(),
        45 * Minute,
        1 * Milli(Metre));
  }

  void FillPositions(Ephemeris<KSP> const& ephemeris,
                     Instant const& initial_time,
                     Time const& duration,
                     std::vector<std::vector<Vector<double, KSP>>>& container) {
    for (Instant t = initial_time;
         t < initial_time + duration;
         t += 45 * Minute) {
      auto const position = [t, &ephemeris](KSPCelestial const& celestial) {
        return ephemeris.trajectory(celestial.body)->
                   EvaluatePosition(t, /*hint=*/nullptr);
      };
      BarycentreCalculator<Position<KSP>, GravitationalParameter>
          jool_system_barycentre;
      for (auto const celestial : jool_system_) {
        jool_system_barycentre.Add(position(*celestial),
                                   celestial->body->gravitational_parameter());
      }
      container.emplace_back();
      for (auto const celestial : jool_system_) {
        container.back().emplace_back(
            (position(*celestial) - jool_system_barycentre.Get()) / Metre);
      }
    }
  };

  Instant const ksp_epoch;

  KSPCelestial sun_;
  KSPCelestial eeloo_;
  KSPCelestial jool_;
  KSPCelestial pol_;
  KSPCelestial bop_;
  KSPCelestial tylo_;
  KSPCelestial vall_;
  KSPCelestial laythe_;
  KSPCelestial dres_;
  KSPCelestial duna_;
  KSPCelestial ike_;
  KSPCelestial kerbin_;
  KSPCelestial minmus_;
  KSPCelestial mun_;
  KSPCelestial eve_;
  KSPCelestial gilly_;
  KSPCelestial moho_;
  std::vector<not_null<KSPCelestial*>> const all_bodies_ = {&sun_,
                                                            &eeloo_,
                                                            &jool_,
                                                            &pol_,
                                                            &bop_,
                                                            &tylo_,
                                                            &vall_,
                                                            &laythe_,
                                                            &dres_,
                                                            &duna_,
                                                            &ike_,
                                                            &kerbin_,
                                                            &minmus_,
                                                            &mun_,
                                                            &eve_,
                                                            &gilly_,
                                                            &moho_};
  std::vector<not_null<KSPCelestial*>> const planets_and_moons_ = {&eeloo_,
                                                                   &jool_,
                                                                   &pol_,
                                                                   &bop_,
                                                                   &tylo_,
                                                                   &vall_,
                                                                   &laythe_,
                                                                   &dres_,
                                                                   &duna_,
                                                                   &ike_,
                                                                   &kerbin_,
                                                                   &minmus_,
                                                                   &mun_,
                                                                   &eve_,
                                                                   &gilly_,
                                                                   &moho_};
  std::vector<not_null<KSPCelestial*>> const jool_system_ =
      {&jool_, &laythe_, &vall_, &tylo_, &pol_, &bop_};
};

TEST_F(KSPSystemTest, KerbalSystem) {
  google::LogToStderr();

  auto const moons = {&laythe_, &vall_, &tylo_, &pol_, &bop_};

  auto const ephemeris = MakeEphemeris();
#if NDEBUG
#if 0
  auto const a_century_hence = ksp_epoch + 100 * JulianYear;
#else  // A small century so the tests don't take too long.
  auto const a_century_hence = ksp_epoch + 5 * JulianYear;
#endif

  LOG(INFO) << "Starting integration";
  ephemeris->Prolong(a_century_hence);
  LOG(INFO) << "Integration done";

  std::map<not_null<KSPCelestial const*>, std::vector<double>>
      extremal_separations_in_m;
  std::map<not_null<KSPCelestial const*>, std::vector<double>> times_in_s;
  Instant t = ksp_epoch;
  std::map<not_null<KSPCelestial const*>, Length> last_separations;
  std::map<not_null<KSPCelestial const*>, Sign> last_separation_changes;

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

  for (auto const* moon : moons) {
    last_separation_changes.emplace(moon, Sign(+1));
  }
  for (int n = 0; t < a_century_hence; ++n, t = ksp_epoch + n * Hour) {
    auto const position = [t, &ephemeris](KSPCelestial const& celestial) {
      return ephemeris->trajectory(celestial.body)->
                 EvaluatePosition(t, /*hint=*/nullptr);
    };
    auto const degrees_of_freedom = [t, &ephemeris](
        KSPCelestial const& celestial) {
      return ephemeris->trajectory(celestial.body)->
                 EvaluateDegreesOfFreedom(t, /*hint=*/nullptr);
    };
    auto const jool_position = position(jool_);

    for (auto const* moon : moons) {
      Length const separation = (jool_position - position(*moon)).Norm();
      Sign const separation_change = Sign(separation - last_separations[moon]);
      if (separation_change != last_separation_changes.at(moon)) {
        extremal_separations_in_m[moon].emplace_back(last_separations[moon] /
                                                     Metre);
        times_in_s[moon].emplace_back((t - 1 * Hour - ksp_epoch) / Second);
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
          KeplerOrbit<KSP>(*jool_.body,
                           MasslessBody(),
                           degrees_of_freedom(bop_) - degrees_of_freedom(jool_),
                           t).elements_at_epoch();
      bop_eccentricities.emplace_back(bop_elements.eccentricity);
      bop_inclinations_in_degrees.emplace_back(bop_elements.inclination /
                                               Degree);
      bop_nodes_in_degrees.emplace_back(
          bop_elements.longitude_of_ascending_node / Degree);
      bop_arguments_of_periapsis_in_degrees.emplace_back(
          bop_elements.argument_of_periapsis / Degree);
    }

    {
      BarycentreCalculator<DegreesOfFreedom<KSP>, GravitationalParameter>
          innermost_jool_system;
      for (auto const* celestial : {&jool_, &laythe_, &vall_, &tylo_}) {
        innermost_jool_system.Add(degrees_of_freedom(*celestial),
                                  celestial->body->gravitational_parameter());
      }
      auto const bop_jacobi_elements =
          KeplerOrbit<KSP>(
              MassiveBody(innermost_jool_system.weight()),
              *bop_.body,
              degrees_of_freedom(bop_) - innermost_jool_system.Get(),
              t).elements_at_epoch();
      bop_jacobi_eccentricities.emplace_back(bop_jacobi_elements.eccentricity);
      bop_jacobi_inclinations_in_degrees.emplace_back(
          bop_jacobi_elements.inclination / Degree);
      bop_jacobi_nodes_in_degrees.emplace_back(
          bop_jacobi_elements.longitude_of_ascending_node / Degree);
      bop_jacobi_arguments_of_periapsis_in_degrees.emplace_back(
          bop_jacobi_elements.argument_of_periapsis / Degree);
    }
  }

  std::vector<std::vector<Vector<double, KSP>>> barycentric_positions_1_year;
  FillPositions(*ephemeris,
                ksp_epoch,
                1 * JulianYear,
                barycentric_positions_1_year);
  std::vector<std::vector<Vector<double, KSP>>> barycentric_positions_2_year;
  FillPositions(*ephemeris,
                ksp_epoch,
                2 * JulianYear,
                barycentric_positions_2_year);

  for (auto const& body_positions : barycentric_positions_1_year) {
    for (auto const& body_position : body_positions) {
      EXPECT_THAT(body_position.Norm(), Lt(3e8));
    }
  }

  std::ofstream file;
  file.open("ksp_system.generated.wl");
  file << mathematica::Assign("laytheTimes", times_in_s[&laythe_]);
  file << mathematica::Assign("vallTimes", times_in_s[&vall_]);
  file << mathematica::Assign("tyloTimes", times_in_s[&tylo_]);
  file << mathematica::Assign("polTimes", times_in_s[&pol_]);
  file << mathematica::Assign("bopTimes", times_in_s[&bop_]);
  file << mathematica::Assign("laytheSeparations",
                              extremal_separations_in_m[&laythe_]);
  file << mathematica::Assign("vallSeparations",
                              extremal_separations_in_m[&vall_]);
  file << mathematica::Assign("tyloSeparations",
                              extremal_separations_in_m[&tylo_]);
  file << mathematica::Assign("polSeparations",
                              extremal_separations_in_m[&pol_]);
  file << mathematica::Assign("bopSeparations",
                              extremal_separations_in_m[&bop_]);

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
  file.close();
#endif
}

}  // namespace physics
}  // namespace principia
