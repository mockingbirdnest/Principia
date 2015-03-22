#pragma once

#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>

#include "glog/logging.h"
#include "quantities/quantities.hpp"

#ifdef ADVANCE_ΔQSTAGE
#error ADVANCE_ΔQSTAGE already defined
#else
#define ADVANCE_ΔQSTAGE(step)                                              \
  do {                                                                     \
    Time const step_evaluated = (step);                                    \
    compute_velocity(p_stage, &v);                                         \
    for (int k = 0; k < dimension; ++k) {                                  \
      Position const Δq = (*Δqstage_previous)[k] + step_evaluated * v[k];  \
      q_stage[k] = q_last[k].value + Δq;                                   \
      (*Δqstage_current)[k] = Δq;                                          \
    }                                                                      \
  } while (false)
#endif

#ifdef ADVANCE_ΔPSTAGE
#error ADVANCE_ΔPSTAGE already defined
#else
#define ADVANCE_ΔPSTAGE(step, q_clock)                                     \
  do {                                                                     \
    Time const step_evaluated = (step);                                    \
    compute_force((q_clock), q_stage, &f);                                 \
    for (int k = 0; k < dimension; ++k) {                                  \
      Momentum const Δp = (*Δpstage_previous)[k] + step_evaluated * f[k];  \
      p_stage[k] = p_last[k].value + Δp;                                   \
      (*Δpstage_current)[k] = Δp;                                          \
    }                                                                      \
  } while (false)
#endif

namespace principia {

using quantities::Quotient;

namespace integrators {

inline SPRKIntegrator const& Leapfrog() {
  static SPRKIntegrator const integrator({ 0.5, 0.5}, { 0.0, 1.0});
  return integrator;
}
inline SPRKIntegrator const& PseudoLeapfrog() {
  static SPRKIntegrator const integrator({ 1.0, 0.0}, { 0.5, 0.5});
  return integrator;
}

inline SPRKIntegrator const& McLachlanAtela1992Order2Optimal() {
  static SPRKIntegrator const integrator({ 0.7071067811865475244,
                                           0.2928932188134524756},
                                         { 0.2928932188134524756,
                                           0.7071067811865475244});
  return integrator;
}

inline SPRKIntegrator const& McLachlan1995S2() {
  static SPRKIntegrator const integrator(
      {0.19318332750378357396, 0.61363334499243285207, 0.19318332750378357396},
      {0.0, 0.5, 0.5});
  return integrator;
}

inline SPRKIntegrator const& Ruth1983() {
  static SPRKIntegrator const integrator({ 2. / 3., -2. / 3., 1.},
                                         { 7. / 24., 3. / 4., -1. / 24.});
  return integrator;
}

inline SPRKIntegrator const& McLachlanAtela1992Order3Optimal() {
  static SPRKIntegrator const integrator({ 0.9196615230173998571,
                                          -0.1879916187991597820,
                                           0.2683300957817599250},
                                         { 0.2683300957817599250,
                                          -0.1879916187991597820,
                                           0.9196615230173998571});
  return integrator;
}

inline SPRKIntegrator const&
CandyRozmus1991ForestRuth1990SynchronousMomenta() {
  static SPRKIntegrator const integrator({ 0.6756035959798288170,
                                          -0.1756035959798288170,
                                          -0.1756035959798288170,
                                           0.6756035959798288170},
                                         { 0.0,
                                           1.351207191959657634,
                                          -1.702414383919315268,
                                           1.351207191959657634});
  return integrator;
}

inline SPRKIntegrator const&
CandyRozmus1991ForestRuth1990SynchronousPositions() {
  static SPRKIntegrator const integrator({ 1.3512071919596576340,
                                          -1.7024143839193152681,
                                           1.3512071919596576340,
                                           0.0},
                                         { 0.6756035959798288170,
                                          -0.1756035959798288170,
                                          -0.1756035959798288170,
                                           0.6756035959798288170});
  return integrator;
}

inline SPRKIntegrator const& McLachlan1995SS5() {
  static SPRKIntegrator const integrator({ 0.14,
                                           0.45273321423383502251,
                                          -0.0927332142338350225,
                                          -0.0927332142338350225,
                                           0.45273321423383502251,
                                           0.14},
                                         { 0.0,
                                           0.28,
                                           0.6254664284676700450,
                                          -0.8109328569353400900,
                                           0.6254664284676700450,
                                           0.28});
  return integrator;
}

inline SPRKIntegrator const& McLachlan1995S4() {
  static SPRKIntegrator const integrator({ 0.16913927992207204518,
                                          -0.29918620390405079951,
                                           1.2600938479639575087,
                                          -0.29918620390405079951,
                                           0.16913927992207204518},
                                         { 0.0,
                                           0.54545454545454545455,
                                          -0.045454545454545454545,
                                          -0.045454545454545454545,
                                           0.54545454545454545455});
  return integrator;
}

inline SPRKIntegrator const& Suzuki1990() {
  static SPRKIntegrator const integrator({ 0.20724538589718786857,
                                           0.41449077179437573714,
                                          -0.12173615769156360571,
                                          -0.12173615769156360571,
                                           0.41449077179437573714,
                                           0.20724538589718786857},
                                         { 0.0,
                                           0.41449077179437573714,
                                           0.41449077179437573714,
                                          -0.65796308717750294857,
                                           0.41449077179437573714,
                                           0.41449077179437573714});
  return integrator;
}

inline SPRKIntegrator const& McLachlan1995S5() {
  static SPRKIntegrator const integrator({ 0.089269454226475244887,
                                          -0.097336042636895508015,
                                           0.50806658841042026313,
                                           0.50806658841042026313,
                                          -0.097336042636895508015,
                                           0.089269454226475244887},
                                         { 0.0,
                                           0.4,
                                          -0.1,
                                           0.4,
                                          -0.1,
                                           0.4});
  return integrator;
}

inline SPRKIntegrator const& BlanesMoan2002S6() {
  static SPRKIntegrator const integrator({ 0.079203696431195700,
                                           0.35317290604977400,
                                          -0.042065080357719500,
                                           0.21937695575349960,
                                          -0.042065080357719500,
                                           0.35317290604977400,
                                           0.079203696431195700},
                                         { 0.0,
                                           0.20951510661336200,
                                          -0.14385177317981800,
                                           0.43433666656645600,
                                           0.43433666656645600,
                                          -0.14385177317981800,
                                           0.20951510661336200});
  return integrator;
}

inline SPRKIntegrator const& Yoshida1990Order6A() {
  static SPRKIntegrator const integrator({ 0.78451361047755726382,
                                           0.23557321335935813369,
                                          -1.17767998417887100695,
                                           1.31518632068391121889,
                                          -1.17767998417887100695,
                                           0.23557321335935813369,
                                           0.78451361047755726382,
                                           0.0},
                                         { 0.392256805238778631910,
                                           0.51004341191845769875,
                                          -0.47105338540975643663,
                                           0.06875316825252010597,
                                           0.06875316825252010597,
                                          -0.47105338540975643663,
                                           0.51004341191845769875,
                                           0.392256805238778631910});
  return integrator;
}

inline SPRKIntegrator const& Yoshida1990Order6B() {
  static SPRKIntegrator const integrator({ 1.43984816797678309102,
                                           0.00426068187079201616,
                                          -2.13228522200145152088,
                                           2.37635274430775282740,
                                          -2.13228522200145152088,
                                           0.00426068187079201616,
                                           1.43984816797678309102,
                                           0.0},
                                         { 0.71992408398839154551,
                                           0.72205442492378755359,
                                          -1.06401227006532975236,
                                           0.12203376115315065326,
                                           0.12203376115315065326,
                                          -1.06401227006532975236,
                                           0.72205442492378755359,
                                           0.71992408398839154551});
  return integrator;
}

inline SPRKIntegrator const& Yoshida1990Order6C() {
  static SPRKIntegrator const integrator({ 1.44778256239929793290,
                                          -2.14403531630538931060,
                                           0.00152886228424927025338,
                                           2.38944778324368421490,
                                           0.00152886228424927025338,
                                          -2.14403531630538931060,
                                           1.44778256239929793290,
                                           0.0},
                                         { 0.72389128119964896645,
                                          -0.34812637695304568885,
                                          -1.07125322701057002017,
                                           1.19548832276396674257,
                                           1.19548832276396674257,
                                          -1.07125322701057002017,
                                          -0.34812637695304568885,
                                           0.72389128119964896645});
  return integrator;
}

inline SPRKIntegrator const& McLachlan1995SS9() {
  static SPRKIntegrator const integrator({ 0.09335,
                                           0.37109851185623919958,
                                           0.34248198631297687861,
                                          -0.35689933723712962525,
                                           0.0499688390679135471,
                                           0.0499688390679135471,
                                          -0.35689933723712962525,
                                           0.34248198631297687861,
                                           0.37109851185623919958,
                                           0.09335},
                                         { 0.0,
                                           0.1867,
                                           0.55549702371247839916,
                                           0.12946694891347535806,
                                          -0.84326562338773460855,
                                           0.9432033015235617027,
                                          -0.84326562338773460855,
                                           0.12946694891347535806,
                                           0.55549702371247839916,
                                           0.1867});
  return integrator;
}

inline SPRKIntegrator const& BlanesMoan2002S10() {
  static SPRKIntegrator const integrator({ 0.050262764400392200,
                                           0.41351430042834400,
                                           0.045079889794397700,
                                          -0.18805485381956900,
                                           0.54196067845078000,
                                          -0.7255255585086898,
                                           0.54196067845078000,
                                          -0.18805485381956900,
                                           0.045079889794397700,
                                           0.41351430042834400,
                                           0.050262764400392200},
                                         { 0.0,
                                           0.14881644790104200,
                                          -0.13238586576778400,
                                           0.067307604692185000,
                                           0.43266640257817500,
                                          -0.01640458940361800,
                                          -0.01640458940361800,
                                           0.43266640257817500,
                                           0.067307604692185000,
                                          -0.13238586576778400,
                                           0.14881644790104200});
  return integrator;
}

inline SPRKIntegrator const& Yoshida1990Order8A() {
  static SPRKIntegrator const integrator({ 1.04242620869970426435,
                                           1.82020630970698006933,
                                           0.157739928123708321343,
                                           2.44002732616634406382,
                                          -0.00716989419709533209981,
                                          -2.44699182370424588929,
                                          -1.61582374150065378479,
                                          -1.7808286265894834253,
                                          -1.61582374150065378479,
                                          -2.44699182370424588929,
                                          -0.00716989419709533209981,
                                           2.44002732616634406382,
                                           0.157739928123708321343,
                                           1.82020630970698006933,
                                           1.04242620869970426435,
                                           0.0},
                                         { 0.521213104349852132174,
                                           1.43131625920334216684,
                                           0.988973118915344195337,
                                           1.29888362714502619258,
                                           1.21642871598462436586,
                                          -1.22708085895067061070,
                                          -2.03140778260244983704,
                                          -1.6983261840450686051,
                                          -1.6983261840450686051,
                                          -2.03140778260244983704,
                                          -1.22708085895067061070,
                                           1.21642871598462436586,
                                           1.29888362714502619258,
                                           0.988973118915344195337,
                                           1.43131625920334216684,
                                           0.521213104349852132174});
  return integrator;
}

inline SPRKIntegrator const& Yoshida1990Order8B() {
  static SPRKIntegrator const integrator({ 1.48819229202921310080,
                                          -2.33864815101041943098,
                                           2.89105148972198900311,
                                          -2.89688250330423987105,
                                           0.00378039588362668223674,
                                           2.89195744315817391244,
                                          -0.00169248587771706559145,
                                          -3.0755169612012526619,
                                          -0.00169248587771706559145,
                                           2.89195744315817391244,
                                           0.00378039588362668223674,
                                          -2.89688250330423987105,
                                           2.89105148972198900311,
                                          -2.33864815101041943098,
                                           1.48819229202921310080,
                                           0.0},
                                         { 0.744096146014606550401,
                                          -0.42522792949060316509,
                                           0.27620166935578478606,
                                          -0.00291550679112543397,
                                          -1.44655105371030659441,
                                           1.44786891952090029734,
                                           1.44513247864022842343,
                                          -1.5386047235394848638,
                                          -1.5386047235394848638,
                                           1.44513247864022842343,
                                           1.44786891952090029734,
                                          -1.44655105371030659441,
                                          -0.00291550679112543397,
                                           0.27620166935578478606,
                                          -0.42522792949060316509,
                                           0.744096146014606550401});
  return integrator;
}

inline SPRKIntegrator const& Yoshida1990Order8C() {
  static SPRKIntegrator const integrator({ 0.629030650210427818049,
                                           1.36934946416874222370,
                                          -1.06458714789183904181,
                                           1.66335809963311356298,
                                          -1.67896928259637402925,
                                          -1.55946803821449795876,
                                           0.311790812418431890510,
                                           1.6589908845439910692,
                                           0.311790812418431890510,
                                          -1.55946803821449795876,
                                          -1.67896928259637402925,
                                           1.66335809963311356298,
                                          -1.06458714789183904181,
                                           1.36934946416874222370,
                                           0.629030650210427818049,
                                           0.0},
                                         { 0.314515325105213909024,
                                           0.999190057189585020872,
                                           0.15238115813845159094,
                                           0.29938547587063726059,
                                          -0.00780559148163023314,
                                          -1.61921866040543599400,
                                          -0.623838612898033034124,
                                           0.98539084848121147984,
                                           0.98539084848121147984,
                                          -0.623838612898033034124,
                                          -1.61921866040543599400,
                                          -0.00780559148163023314,
                                           0.29938547587063726059,
                                           0.15238115813845159094,
                                           0.999190057189585020872,
                                           0.314515325105213909024});
  return integrator;
}

inline SPRKIntegrator const& Yoshida1990Order8D() {
  static SPRKIntegrator const integrator({ 0.914844246229642658287,
                                           0.253693336566286009974,
                                          -1.44485223686030647660,
                                          -0.158240635368502468458,
                                           1.93813913762291232471,
                                          -1.96061023297558163691,
                                           0.102799849392219431139,
                                           1.7084530707866603157,
                                           0.102799849392219431139,
                                          -1.96061023297558163691,
                                           1.93813913762291232471,
                                          -0.158240635368502468458,
                                          -1.44485223686030647660,
                                           0.253693336566286009974,
                                           0.914844246229642658287,
                                           0.0},
                                         { 0.457422123114821329143,
                                           0.584268791397964334130,
                                          -0.595579450147010233314,
                                          -0.801546436114404472530,
                                           0.88994925112720492813,
                                          -0.01123554767633465610,
                                          -0.92890519179168110289,
                                           0.90562646008943987343,
                                           0.90562646008943987343,
                                          -0.92890519179168110289,
                                          -0.01123554767633465610,
                                           0.88994925112720492813,
                                          -0.801546436114404472530,
                                          -0.595579450147010233314,
                                           0.584268791397964334130,
                                           0.457422123114821329143});
  return integrator;
}

inline SPRKIntegrator const& Yoshida1990Order8E() {
  static SPRKIntegrator const integrator({ 1.30300165757516838484,
                                           0.107990467718098279648,
                                          -2.04809795883490205633,
                                           0.00536018921375238082832,
                                          -0.0719180053650705075005,
                                           2.52778927318028339169,
                                           0.0227738840126312259937,
                                          -2.6937990149999221983,
                                           0.0227738840126312259937,
                                           2.52778927318028339169,
                                          -0.0719180053650705075005,
                                           0.00536018921375238082832,
                                          -2.04809795883490205633,
                                           0.107990467718098279648,
                                           1.30300165757516838484,
                                           0.0},
                                         { 0.651500828787584192418,
                                           0.705496062646633332241,
                                          -0.97005374555840188834,
                                          -1.02136888481057483775,
                                          -0.0332789080756590633361,
                                           1.22793563390760644210,
                                           1.27528157859645730884,
                                          -1.33551256549364548617,
                                          -1.33551256549364548617,
                                           1.27528157859645730884,
                                           1.22793563390760644210,
                                          -0.0332789080756590633361,
                                          -1.02136888481057483775,
                                          -0.97005374555840188834,
                                           0.705496062646633332241,
                                           0.651500828787584192418});
  return integrator;
}

inline SPRKIntegrator const& McLachlan1995SS15() {
  static SPRKIntegrator const integrator({ 0.37083518217530647672,
                                           0.1662847692752906797,
                                          -0.1091730577518966070,
                                          -0.1915538804099219434,
                                          -0.13739914490621317141,
                                           0.3168445497744770538,
                                           0.32495900532103239021,
                                          -0.2407974234780748787,
                                          -0.2407974234780748787,
                                           0.32495900532103239021,
                                           0.3168445497744770538,
                                          -0.13739914490621317141,
                                          -0.1915538804099219434,
                                          -0.1091730577518966070,
                                           0.1662847692752906797,
                                           0.37083518217530647672},
                                         { 0.0,
                                           0.7416703643506129534,
                                          -0.4091008258000315940,
                                           0.1907547102962383800,
                                          -0.5738624711160822667,
                                           0.2990641813036559238,
                                           0.3346249182452981838,
                                           0.3152930923967665966,
                                          -0.796887939352916354,
                                           0.3152930923967665966,
                                           0.3346249182452981838,
                                           0.2990641813036559238,
                                          -0.5738624711160822667,
                                           0.1907547102962383800,
                                          -0.4091008258000315940,
                                           0.7416703643506129534});
  return integrator;
}

inline SPRKIntegrator const& McLachlan1995SS17() {
  static SPRKIntegrator const integrator({ 0.064432989690721649485,
                                           0.35519003324334713070,
                                           0.0856693578177004124,
                                          -0.11251421787663120244,
                                          -0.1122027038521318433,
                                          -0.13257320117041968914,
                                           0.2113707207368458462,
                                           0.2966460921549872546,
                                          -0.1560190707444195585,
                                          -0.1560190707444195585,
                                           0.2966460921549872546,
                                           0.2113707207368458462,
                                          -0.13257320117041968914,
                                          -0.1122027038521318433,
                                          -0.11251421787663120244,
                                           0.0856693578177004124,
                                           0.35519003324334713070,
                                           0.064432989690721649485},
                                         { 0.0,
                                           0.12886597938144329897,
                                           0.5815140871052509624,
                                          -0.41017537146985013753,
                                           0.1851469357165877327,
                                          -0.40955234342085141934,
                                           0.14440594108001204106,
                                           0.27833550039367965131,
                                           0.31495668391629485789,
                                          -0.626994825405133975,
                                           0.31495668391629485789,
                                           0.27833550039367965131,
                                           0.14440594108001204106,
                                          -0.40955234342085141934,
                                           0.1851469357165877327,
                                          -0.41017537146985013753,
                                           0.5815140871052509624,
                                           0.12886597938144329897});
  return integrator;
}

inline SPRKIntegrator::SPRKIntegrator(std::vector<double> const& a,
                                      std::vector<double> const& b)
    : SRKNIntegrator(a, b) {}

template<typename Position, typename Momentum,
         typename AutonomousRightHandSideComputation,
         typename RightHandSideComputation>
void SPRKIntegrator::SolveIncrement(
    RightHandSideComputation compute_force,
    AutonomousRightHandSideComputation compute_velocity,
    Parameters<Position, Momentum> const& parameters,
    not_null<Solution<Position, Momentum>*> const solution) const {
  switch (vanishing_coefficients_) {
    case kNone:
      SolveIncrementOptimized<kNone>(
          compute_force, compute_velocity, parameters, solution);
      break;
    case kFirstBVanishes:
      SolveIncrementOptimized<kFirstBVanishes>(
          compute_force, compute_velocity, parameters, solution);
      break;
    case kLastAVanishes:
      SolveIncrementOptimized<kLastAVanishes>(
          compute_force, compute_velocity, parameters, solution);
      break;
    default:
      LOG(FATAL) << "Invalid vanishing coefficients";
  }
}

template<SPRKIntegrator::VanishingCoefficients vanishing_coefficients,
         typename Position, typename Momentum,
         typename AutonomousRightHandSideComputation,
         typename RightHandSideComputation>
void SPRKIntegrator::SolveIncrementOptimized(
    RightHandSideComputation compute_force,
    AutonomousRightHandSideComputation compute_velocity,
    Parameters<Position, Momentum> const& parameters,
    not_null<Solution<Position, Momentum>*> const solution) const {
  int const dimension = parameters.initial.positions.size();

  std::vector<Position> Δqstage0(dimension);
  std::vector<Position> Δqstage1(dimension);
  std::vector<Momentum> Δpstage0(dimension);
  std::vector<Momentum> Δpstage1(dimension);
  std::vector<Position>* Δqstage_current = &Δqstage1;
  std::vector<Position>* Δqstage_previous = &Δqstage0;
  std::vector<Momentum>* Δpstage_current = &Δpstage1;
  std::vector<Momentum>* Δpstage_previous = &Δpstage0;

  // Dimension the result.
  int const capacity = parameters.sampling_period == 0 ?
    1 :
    static_cast<int>(
        ceil((((parameters.tmax - parameters.initial.time.value) /
                    parameters.Δt) + 1) /
                parameters.sampling_period)) + 1;
  solution->clear();
  solution->reserve(capacity);

  std::vector<DoublePrecision<Position>> q_last(parameters.initial.positions);
  std::vector<DoublePrecision<Momentum>> p_last(parameters.initial.momenta);
  int sampling_phase = 0;

  std::vector<Position> q_stage(dimension);
  std::vector<Momentum> p_stage(dimension);
  std::vector<Quotient<Momentum, Time>> f(dimension);  // Current forces.
  std::vector<Quotient<Position, Time>> v(dimension);  // Current velocities.

  // The following quantity is generally equal to |Δt|, but during the last
  // iteration, if |tmax_is_exact|, it may differ significantly from |Δt|.
  Time h = parameters.Δt;  // Constant for now.

  // During one iteration of the outer loop below we process the time interval
  // [|tn|, |tn| + |h|[.  |tn| is computed using compensated summation to make
  // sure that we don't have drifts.
  DoublePrecision<Time> tn = parameters.initial.time;

  // Whether position and momentum are synchronized between steps, relevant for
  // first-same-as-last (FSAL) integrators. Time is always synchronous with
  // position.
  bool q_and_p_are_synchronized = true;
  bool should_synchronize = false;

  // Integration.  For details see Wolfram Reference,
  // http://reference.wolfram.com/mathematica/tutorial/NDSolveSPRK.html#74387056
  bool at_end = !parameters.tmax_is_exact && parameters.tmax < tn.value + h;
  while (!at_end) {
    // Check if this is the last interval and if so process it appropriately.
    if (parameters.tmax_is_exact) {
      // If |tn| is getting close to |tmax|, use |tmax| as the upper bound of
      // the interval and update |h| accordingly.  The bound chosen here for
      // |tmax| ensures that we don't end up with a ridiculously small last
      // interval: we'd rather make the last interval a bit bigger.  More
      // precisely, the last interval generally has a length between 0.5 Δt and
      // 1.5 Δt, unless it is also the first interval.
      // NOTE(phl): This may lead to convergence as bad as (1.5 Δt)^5 rather
      // than Δt^5.
      if (parameters.tmax <= tn.value + 3 * h / 2) {
        at_end = true;
        h = (parameters.tmax - tn.value) - tn.error;
      }
    } else if (parameters.tmax < tn.value + 2 * h) {
      // If the next interval would overshoot, make this the last interval but
      // stick to the same step.
      at_end = true;
    }
    // Here |h| is the length of the current time interval and |tn| is its
    // start.

    // Increment SPRK step from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 3.
    for (int k = 0; k < dimension; ++k) {
      (*Δqstage_current)[k] = Position();
      (*Δpstage_current)[k] = Momentum();
      q_stage[k] = q_last[k].value;
    }

    if (vanishing_coefficients != kNone) {
      should_synchronize = at_end ||
                           (parameters.sampling_period != 0 &&
                            sampling_phase % parameters.sampling_period == 0);
    }

    if (vanishing_coefficients == kFirstBVanishes &&
        q_and_p_are_synchronized) {
      // Desynchronize.
      std::swap(Δqstage_current, Δqstage_previous);
      for (int k = 0; k < dimension; ++k) {
        p_stage[k] = p_last[k].value;
      }
      ADVANCE_ΔQSTAGE(first_same_as_last_->first * h);
      q_and_p_are_synchronized = false;
    }
    for (int i = 0; i < stages_; ++i) {
      std::swap(Δqstage_current, Δqstage_previous);
      std::swap(Δpstage_current, Δpstage_previous);

      // Beware, the p/q order matters here, the two computations depend on one
      // another.

      // By using |tn.error| below we get a time value which is possibly a wee
      // bit more precise.
      if (vanishing_coefficients == kLastAVanishes &&
          q_and_p_are_synchronized && i == 0) {
        ADVANCE_ΔPSTAGE(first_same_as_last_->first * h,
                        tn.value);
        q_and_p_are_synchronized = false;
      } else {
        ADVANCE_ΔPSTAGE(b_[i] * h, tn.value + (tn.error + c_[i] * h));
      }

      if (vanishing_coefficients == kFirstBVanishes &&
          should_synchronize && i == stages_ - 1) {
        ADVANCE_ΔQSTAGE(first_same_as_last_->last * h);
        q_and_p_are_synchronized = true;
      } else {
        ADVANCE_ΔQSTAGE(a_[i] * h);
      }
    }
    if (vanishing_coefficients == kLastAVanishes && should_synchronize) {
      std::swap(Δpstage_current, Δpstage_previous);
      // TODO(egg): the second parameter below is really just tn.value + h.
      ADVANCE_ΔPSTAGE(first_same_as_last_->last * h,
                      tn.value + h);
      q_and_p_are_synchronized = true;
    }
    // Compensated summation from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 2.
    for (int k = 0; k < dimension; ++k) {
      q_last[k].Increment((*Δqstage_current)[k]);
      p_last[k].Increment((*Δpstage_current)[k]);
      q_stage[k] = q_last[k].value;
      p_stage[k] = p_last[k].value;
    }
    tn.Increment(h);

    if (parameters.sampling_period != 0) {
      if (sampling_phase % parameters.sampling_period == 0) {
        solution->emplace_back();
        SystemState<Position, Momentum>* state = &solution->back();
        state->time = tn;
        state->positions.reserve(dimension);
        state->momenta.reserve(dimension);
        for (int k = 0; k < dimension; ++k) {
          state->positions.emplace_back(q_last[k]);
          state->momenta.emplace_back(p_last[k]);
        }
      }
      ++sampling_phase;
    }
  }

  if (parameters.sampling_period == 0) {
    solution->emplace_back();
    SystemState<Position, Momentum>* state = &solution->back();
    state->time = tn;
    state->positions.reserve(dimension);
    state->momenta.reserve(dimension);
    for (int k = 0; k < dimension; ++k) {
      state->positions.emplace_back(q_last[k]);
      state->momenta.emplace_back(p_last[k]);
    }
  }
}

}  // namespace integrators
}  // namespace principia

#undef ADVANCE_ΔQSTAGE
#undef ADVANCE_ΔPSTAGE
