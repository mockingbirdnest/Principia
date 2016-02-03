
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
#include "testing_utilities/numerics.hpp"

namespace principia {

using quantities::astronomy::JulianYear;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using testing_utilities::RelativeError;

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
    sun.owned_body = make_not_null_unique<MassiveBody>(
        1.1723327948324908E+18 * SIUnit<GravitationalParameter>());
    eeloo.owned_body = make_not_null_unique<MassiveBody>(
        74410814527.049576 * SIUnit<GravitationalParameter>());
    jool.owned_body = make_not_null_unique<MassiveBody>(
        282528004209995.31 * SIUnit<GravitationalParameter>());
    pol.owned_body = make_not_null_unique<MassiveBody>(
        721702080.00000012 * SIUnit<GravitationalParameter>());
    bop.owned_body = make_not_null_unique<MassiveBody>(
        2486834944.414907 * SIUnit<GravitationalParameter>());
    tylo.owned_body = make_not_null_unique<MassiveBody>(
        2825280042099.9531 * SIUnit<GravitationalParameter>());
    vall.owned_body = make_not_null_unique<MassiveBody>(
        207481499473.75098 * SIUnit<GravitationalParameter>());
    laythe.owned_body = make_not_null_unique<MassiveBody>(
        1962000029236.0784 * SIUnit<GravitationalParameter>());
    dres.owned_body = make_not_null_unique<MassiveBody>(
        21484488600.000004 * SIUnit<GravitationalParameter>());
    duna.owned_body = make_not_null_unique<MassiveBody>(
        301363211975.09772 * SIUnit<GravitationalParameter>());
    ike.owned_body = make_not_null_unique<MassiveBody>(
        18568368573.144012 * SIUnit<GravitationalParameter>());
    kerbin.owned_body = make_not_null_unique<MassiveBody>(
        3531600000000 * SIUnit<GravitationalParameter>());
    minmus.owned_body = make_not_null_unique<MassiveBody>(
        1765800026.3124719 * SIUnit<GravitationalParameter>());
    mun.owned_body = make_not_null_unique<MassiveBody>(
        65138397520.780701 * SIUnit<GravitationalParameter>());
    eve.owned_body = make_not_null_unique<MassiveBody>(
        8171730229210.874 * SIUnit<GravitationalParameter>());
    gilly.owned_body = make_not_null_unique<MassiveBody>(
        8289449.814716354 * SIUnit<GravitationalParameter>());
    moho.owned_body = make_not_null_unique<MassiveBody>(
        168609378654.50949 * SIUnit<GravitationalParameter>());
    for (auto const celestial : all_bodies) {
      celestial->body = celestial->owned_body.get();
    }
    eeloo.parent = &sun;
    jool.parent = &sun;
    pol.parent = &jool;
    bop.parent = &jool;
    tylo.parent = &jool;
    vall.parent = &jool;
    laythe.parent = &jool;
    dres.parent = &sun;
    duna.parent = &sun;
    ike.parent = &duna;
    kerbin.parent = &sun;
    minmus.parent = &kerbin;
    mun.parent = &kerbin;
    eve.parent = &sun;
    gilly.parent = &eve;
    moho.parent = &sun;
    eeloo.elements.eccentricity = +2.60000000000000009e-01;
    eeloo.elements.mean_motion = +4.00223155970064009e-08 * (Radian / Second);
    eeloo.elements.inclination = +1.07337748997651278e-01 * Radian;
    eeloo.elements.longitude_of_ascending_node =
        +8.72664625997164767e-01 * Radian;
    eeloo.elements.argument_of_periapsis = +4.53785605518525692e+00 * Radian;
    eeloo.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    jool.elements.eccentricity = +5.00000007450581013e-02;
    jool.elements.mean_motion = +6.00334352457231946e-08 * (Radian / Second);
    jool.elements.inclination = +2.27590937955459392e-02 * Radian;
    jool.elements.longitude_of_ascending_node =
        +9.07571211037051406e-01 * Radian;
    jool.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    jool.elements.mean_anomaly = +1.00000001490115994e-01 * Radian;
    pol.elements.eccentricity = +1.70850000000000002e-01;
    pol.elements.mean_motion = +6.96658945572122982e-06 * (Radian / Second);
    pol.elements.inclination = +7.41764932097590118e-02 * Radian;
    pol.elements.longitude_of_ascending_node =
        +3.49065850398865909e-02 * Radian;
    pol.elements.argument_of_periapsis = +2.61799387799149408e-01 * Radian;
    pol.elements.mean_anomaly = +8.99999976158141979e-01 * Radian;
    bop.elements.eccentricity = +2.34999999403953996e-01;
    bop.elements.mean_motion = +4.64439297048081960e-06 * (Radian / Second);
    bop.elements.inclination = +2.61799387799149408e-01 * Radian;
    bop.elements.longitude_of_ascending_node =
        +1.74532925199432948e-01 * Radian;
    bop.elements.argument_of_periapsis = +4.36332312998582383e-01 * Radian;
    bop.elements.mean_anomaly = +8.99999976158141979e-01 * Radian;
    tylo.elements.eccentricity = +0.00000000000000000e+00;
    tylo.elements.mean_motion = +1.94051054171045988e-05 * (Radian / Second);
    tylo.elements.inclination = +4.36332319500439990e-04 * Radian;
    tylo.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    tylo.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    tylo.elements.mean_anomaly = +3.14159265358979001e+00 * Radian;
    vall.elements.eccentricity = +0.00000000000000000e+00;
    vall.elements.mean_motion = +4.79720588121814983e-05 * (Radian / Second);
    vall.elements.inclination = +0.00000000000000000e+00 * Radian;
    vall.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    vall.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    vall.elements.mean_anomaly = +0.00000000000000000e+00 * Radian;
    laythe.elements.eccentricity = +0.00000000000000000e+00;
    laythe.elements.mean_motion = +1.18593451424947995e-04 * (Radian / Second);
    laythe.elements.inclination = +0.00000000000000000e+00 * Radian;
    laythe.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    laythe.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    laythe.elements.mean_anomaly = +3.14159265358979001e+00 * Radian;
    dres.elements.eccentricity = +1.44999999999999990e-01;
    dres.elements.mean_motion = +1.31191970097993002e-07 * (Radian / Second);
    dres.elements.inclination = +8.72664625997164739e-02 * Radian;
    dres.elements.longitude_of_ascending_node =
        +4.88692190558412243e+00 * Radian;
    dres.elements.argument_of_periapsis = +1.57079632679489656e+00 * Radian;
    dres.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    duna.elements.eccentricity = +5.09999990463256975e-02;
    duna.elements.mean_motion = +3.62866884706430976e-07 * (Radian / Second);
    duna.elements.inclination = +1.04719752778990849e-03 * Radian;
    duna.elements.longitude_of_ascending_node =
        +2.36492113645231639e+00 * Radian;
    duna.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    duna.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    ike.elements.eccentricity = +2.99999993294477012e-02;
    ike.elements.mean_motion = +9.59003407994517016e-05 * (Radian / Second);
    ike.elements.inclination = +3.49065855600351992e-03 * Radian;
    ike.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    ike.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    ike.elements.mean_anomaly = +1.70000004768372004e+00 * Radian;
    kerbin.elements.eccentricity = +0.00000000000000000e+00;
    kerbin.elements.mean_motion = +6.82691894080843017e-07 * (Radian / Second);
    kerbin.elements.inclination = +0.00000000000000000e+00 * Radian;
    kerbin.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    kerbin.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    kerbin.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    minmus.elements.eccentricity = +0.00000000000000000e+00;
    minmus.elements.mean_motion = +5.83228807719951003e-06 * (Radian / Second);
    minmus.elements.inclination = +1.04719755119659780e-01 * Radian;
    minmus.elements.longitude_of_ascending_node =
        +1.36135681655557694e+00 * Radian;
    minmus.elements.argument_of_periapsis = +6.63225115757845263e-01 * Radian;
    minmus.elements.mean_anomaly = +8.99999976158141979e-01 * Radian;
    mun.elements.eccentricity = +0.00000000000000000e+00;
    mun.elements.mean_motion = +4.52078533000627999e-05 * (Radian / Second);
    mun.elements.inclination = +0.00000000000000000e+00 * Radian;
    mun.elements.longitude_of_ascending_node =
        +0.00000000000000000e+00 * Radian;
    mun.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    mun.elements.mean_anomaly = +1.70000004768372004e+00 * Radian;
    eve.elements.eccentricity = +9.99999977648258036e-03;
    eve.elements.mean_motion = +1.11049676511037010e-06 * (Radian / Second);
    eve.elements.inclination = +3.66519126274052684e-02 * Radian;
    eve.elements.longitude_of_ascending_node =
        +2.61799387799149408e-01 * Radian;
    eve.elements.argument_of_periapsis = +0.00000000000000000e+00 * Radian;
    eve.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
    gilly.elements.eccentricity = +5.50000011920928955e-01;
    gilly.elements.mean_motion = +1.61692985452753988e-05 * (Radian / Second);
    gilly.elements.inclination = +2.09439510239319560e-01 * Radian;
    gilly.elements.longitude_of_ascending_node =
        +1.39626340159546358e+00 * Radian;
    gilly.elements.argument_of_periapsis = +1.74532925199432948e-01 * Radian;
    gilly.elements.mean_anomaly = +8.99999976158141979e-01 * Radian;
    moho.elements.eccentricity = +2.00000002980231989e-01;
    moho.elements.mean_motion = +2.83568694188237007e-06 * (Radian / Second);
    moho.elements.inclination = +1.22173047639603072e-01 * Radian;
    moho.elements.longitude_of_ascending_node =
        +1.22173047639603061e+00 * Radian;
    moho.elements.argument_of_periapsis = +2.61799387799149408e-01 * Radian;
    moho.elements.mean_anomaly = +3.14000010490416992e+00 * Radian;
  }

  not_null<std::unique_ptr<Ephemeris<KSP>>> MakeEphemeris() {
    HierarchicalSystem<KSP> hierarchical_system(std::move(sun.owned_body));
    for (auto const celestial : planets_and_moons ) {
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
        integrators::BlanesMoan2002SRKN14A<Position<KSP>>(),
        45 * Minute,
        1 * Milli(Metre));
  }

  Instant const ksp_epoch;

  KSPCelestial sun;
  KSPCelestial eeloo;
  KSPCelestial jool;
  KSPCelestial pol;
  KSPCelestial bop;
  KSPCelestial tylo;
  KSPCelestial vall;
  KSPCelestial laythe;
  KSPCelestial dres;
  KSPCelestial duna;
  KSPCelestial ike;
  KSPCelestial kerbin;
  KSPCelestial minmus;
  KSPCelestial mun;
  KSPCelestial eve;
  KSPCelestial gilly;
  KSPCelestial moho;
  std::vector<not_null<KSPCelestial*>> const all_bodies = {&sun,
                                                           &eeloo,
                                                           &jool,
                                                           &pol,
                                                           &bop,
                                                           &tylo,
                                                           &vall,
                                                           &laythe,
                                                           &dres,
                                                           &duna,
                                                           &ike,
                                                           &kerbin,
                                                           &minmus,
                                                           &mun,
                                                           &eve,
                                                           &gilly,
                                                           &moho};
  std::vector<not_null<KSPCelestial*>> const planets_and_moons = {&eeloo,
                                                                  &jool,
                                                                  &pol,
                                                                  &bop,
                                                                  &tylo,
                                                                  &vall,
                                                                  &laythe,
                                                                  &dres,
                                                                  &duna,
                                                                  &ike,
                                                                  &kerbin,
                                                                  &minmus,
                                                                  &mun,
                                                                  &eve,
                                                                  &gilly,
                                                                  &moho};
};

TEST_F(KSPSystemTest, KerbalSystem) {
  google::LogToStderr();

  bop.elements.inclination = π * Radian - bop.elements.inclination;
  bop.elements.mean_motion = *pol.elements.mean_motion / 0.7;
  bop.elements.eccentricity *= 1.2;
  LOG(INFO) << bop.elements;
  LOG(INFO) << pol.elements;

  auto const moons = {&laythe, &vall, &tylo, &pol, &bop};

  auto const ephemeris = MakeEphemeris();
  auto const a_century_hence = ksp_epoch + 100 * JulianYear;

  LOG(INFO) << "Starting integration";
  ephemeris->Prolong(a_century_hence);
  LOG(INFO) << "Done";

  auto const jool_trajectory = ephemeris->trajectory(jool.body);
  std::map<not_null<KSPCelestial const*>, ContinuousTrajectory<KSP> const*>
      moon_trajectories;
  for (auto const* moon : moons) {
    moon_trajectories[moon] = ephemeris->trajectory(moon->body);
  }
  std::map<not_null<KSPCelestial const*>, std::vector<double>>
      extremal_separations_in_m;
  std::map<not_null<KSPCelestial const*>, std::vector<double>> times_in_s;
  Instant t = ksp_epoch;
  std::map<not_null<KSPCelestial const*>, Length> last_separations;
  std::map<not_null<KSPCelestial const*>, Sign> last_separation_changes;
  for (auto const* moon : moons) {
    last_separation_changes.emplace(moon, Sign(+1));
  }
  for (int n = 0; t < a_century_hence; ++n, t = ksp_epoch + n * Hour) {
    auto const jool_position =
        jool_trajectory->EvaluatePosition(t, /*hint=*/nullptr);
    for (auto const* moon : moons) {
      auto const moon_position =
          moon_trajectories[moon]->EvaluatePosition(t, /*hint=*/nullptr);
      Length const separation = (jool_position - moon_position).Norm();
      Sign separation_change = Sign(separation - last_separations[moon]);
      if (separation_change != last_separation_changes.at(moon)) {
        extremal_separations_in_m[moon].emplace_back(last_separations[moon] /
                                                     Metre);
        times_in_s[moon].emplace_back((t - 1 * Hour - ksp_epoch) / Second);
      }
      last_separations[moon] = separation;
      last_separation_changes.at(moon) = separation_change;
    }
  }
  std::ofstream file;
  file.open("ksp_system.generated.wl");
  file << mathematica::Assign("laytheTimes", times_in_s[&laythe]);
  file << mathematica::Assign("vallTimes", times_in_s[&vall]);
  file << mathematica::Assign("tyloTimes", times_in_s[&tylo]);
  file << mathematica::Assign("polTimes", times_in_s[&pol]);
  file << mathematica::Assign("bopTimes", times_in_s[&bop]);
  file << mathematica::Assign("laytheSeparations",
                              extremal_separations_in_m[&laythe]);
  file << mathematica::Assign("vallSeparations",
                              extremal_separations_in_m[&vall]);
  file << mathematica::Assign("tyloSeparations",
                              extremal_separations_in_m[&tylo]);
  file << mathematica::Assign("polSeparations",
                              extremal_separations_in_m[&pol]);
  file << mathematica::Assign("bopSeparations",
                              extremal_separations_in_m[&bop]);
  file.close();
}

}  // namespace physics
}  // namespace principia
