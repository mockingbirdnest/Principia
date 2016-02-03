
#include <map>
#include <string>
#include <vector>

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

using integrators::McLachlanAtela1992Order5Optimal;
using quantities::astronomy::JulianYear;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
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
    not_null<std::unique_ptr<MassiveBody>> body;
  };

  KSPSystemTest() {
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
};
