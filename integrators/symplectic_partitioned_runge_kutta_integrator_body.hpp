
#pragma once

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include "base/mod.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_runge_kutta_nyström_integrator {

using base::mod;

template<typename Method, typename Position>
SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
    SymplecticPartitionedRungeKuttaIntegrator() {
  // TODO(phl): This might be turned into a static_assert.
  if (Method::first_same_as_last) {
    CHECK_EQ(0.0, Method::a[Method::stages - 1]);
  }
  if (Method::time_reversible) {
    CHECK(Method::first_same_as_last);
    for (int i = 0; i < Method::stages - 1; ++i) {
      CHECK_EQ(Method::a[i], Method::a[Method::stages - 2 - i]);
    }
    for (int i = 0; i < Method::stages; ++i) {
      CHECK_EQ(Method::b[i], Method::b[stages - 1 - i]);
    }
  }
}

template<typename Position,
         typename Momentum,
         int order_,
         bool time_reversible_,
         int evaluations_,
         bool first_same_as_last_>
template<CompositionMethod composition_method>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      order_,
                                      time_reversible_,
                                      evaluations_,
                                      composition_method> const&
    SymplecticPartitionedRungeKuttaIntegrator<
        Position,
        Momentum,
        order_,
        time_reversible_,
        evaluations_,
        first_same_as_last_>::AsRungeKuttaNyströmIntegrator() const {
  using SRKN = SymplecticRungeKuttaNyströmIntegrator<Position,
                                                     order_,
                                                     time_reversible_,
                                                     evaluations_,
                                                     composition_method>;
  static_assert(first_same_as_last
                    ? composition_method == ABA || composition_method == BAB
                    : composition_method == BA,
                "requested |composition_method| inconsistent with the "
                "properties of this integrator");
  // The |reinterpret_cast|s are ugly, but everything else I can think of is a
  // rabbit hole of metaprogramming.  They are not UB: the option that will be
  // picked has the right type.
  std::unique_ptr<SRKN>& method =
      composition_method == BA
          ? reinterpret_cast<std::unique_ptr<SRKN>&>(ba_srkn_)
          : composition_method == ABA
                ? reinterpret_cast<std::unique_ptr<SRKN>&>(aba_srkn_)
                : reinterpret_cast<std::unique_ptr<SRKN>&>(bab_srkn_);
  if (method == nullptr) {
    if (composition_method == ABA) {
      FixedVector<double, stages_> shifted_a;
      // |*this| is a |BAB| method, with A and B interchangeable.  Exchanging A
      // and B shifts |a_| (because |ABA| means b₀ vanishes, whereas |BAB| means
      // aᵣ vanishes).
      for (int i = 0; i < stages_; ++i) {
        shifted_a[i] = a_[mod(i - 1, stages_)];
      }
      method = std::make_unique<SRKN>(
          serialization::FixedStepSizeIntegrator::DUMMY, b_, shifted_a);
    } else {
      method = std::make_unique<SRKN>(
          serialization::FixedStepSizeIntegrator::DUMMY, a_, b_);
    }
  }
  return *method;
}

}  // namespace internal_symplectic_runge_kutta_nyström_integrator

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/2,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/1,
                                          /*first_same_as_last=*/true> const&
NewtonDelambreStørmerVerletLeapfrog() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/2,
      /*time_reversible=*/true,
      /*evaluations=*/1,
      /*first_same_as_last=*/true> const integrator({1.0, 0.0}, {0.5, 0.5});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/3,
                                          /*time_reversible=*/false,
                                          /*evaluations=*/3,
                                          /*first_same_as_last=*/false> const&
Ruth1983() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/3,
      /*time_reversible=*/false,
      /*evaluations=*/3,
      /*first_same_as_last=*/false> const integrator({2. / 3., -2. / 3., 1.},
                                                     {7. / 24.,
                                                      3. / 4.,
                                                      -1. / 24.});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/5,
                                          /*first_same_as_last=*/true> const&
Suzuki1990() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/4,
      /*time_reversible=*/true,
      /*evaluations=*/5,
      /*first_same_as_last=*/true> const integrator({+0.41449077179437573714,
                                                     +0.41449077179437573714,
                                                     -0.65796308717750294857,
                                                     +0.41449077179437573714,
                                                     +0.41449077179437573714,
                                                     +0.0},
                                                    {+0.20724538589718786857,
                                                     +0.41449077179437573714,
                                                     -0.12173615769156360571,
                                                     -0.12173615769156360571,
                                                     +0.41449077179437573714,
                                                     +0.20724538589718786857});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/7,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order6A() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/6,
      /*time_reversible=*/true,
      /*evaluations=*/7,
      /*first_same_as_last=*/true> const integrator({+0.78451361047755726382,
                                                     +0.23557321335935813369,
                                                     -1.17767998417887100695,
                                                     +1.31518632068391121889,
                                                     -1.17767998417887100695,
                                                     +0.23557321335935813369,
                                                     +0.78451361047755726382,
                                                     +0.0},
                                                    {+0.392256805238778631910,
                                                     +0.51004341191845769875,
                                                     -0.47105338540975643663,
                                                     +0.06875316825252010597,
                                                     +0.06875316825252010597,
                                                     -0.47105338540975643663,
                                                     +0.51004341191845769875,
                                                     +0.392256805238778631910});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/7,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order6B() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/6,
      /*time_reversible=*/true,
      /*evaluations=*/7,
      /*first_same_as_last=*/true> const integrator({+1.43984816797678309102,
                                                     +0.00426068187079201616,
                                                     -2.13228522200145152088,
                                                     +2.37635274430775282740,
                                                     -2.13228522200145152088,
                                                     +0.00426068187079201616,
                                                     +1.43984816797678309102,
                                                     +0.0},
                                                    {+0.71992408398839154551,
                                                     +0.72205442492378755359,
                                                     -1.06401227006532975236,
                                                     +0.12203376115315065326,
                                                     +0.12203376115315065326,
                                                     -1.06401227006532975236,
                                                     +0.72205442492378755359,
                                                     +0.71992408398839154551});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/7,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order6C() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/6,
      /*time_reversible=*/true,
      /*evaluations=*/7,
      /*first_same_as_last=*/true> const integrator({+1.44778256239929793290,
                                                     -2.14403531630538931060,
                                                     +0.00152886228424927025338,
                                                     +2.38944778324368421490,
                                                     +0.00152886228424927025338,
                                                     -2.14403531630538931060,
                                                     +1.44778256239929793290,
                                                     +0.0},
                                                    {+0.72389128119964896645,
                                                     -0.34812637695304568885,
                                                     -1.07125322701057002017,
                                                     +1.19548832276396674257,
                                                     +1.19548832276396674257,
                                                     -1.07125322701057002017,
                                                     -0.34812637695304568885,
                                                     +0.72389128119964896645});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8A() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/8,
      /*time_reversible=*/true,
      /*evaluations=*/15,
      /*first_same_as_last=*/true> const integrator({+1.04242620869970426435,
                                                     +1.82020630970698006933,
                                                     +0.157739928123708321343,
                                                     +2.44002732616634406382,
                                                     -0.00716989419709533209981,
                                                     -2.44699182370424588929,
                                                     -1.61582374150065378479,
                                                     -1.7808286265894834253,
                                                     -1.61582374150065378479,
                                                     -2.44699182370424588929,
                                                     -0.00716989419709533209981,
                                                     +2.44002732616634406382,
                                                     +0.157739928123708321343,
                                                     +1.82020630970698006933,
                                                     +1.04242620869970426435,
                                                     +0.0},
                                                    {+0.521213104349852132174,
                                                     +1.43131625920334216684,
                                                     +0.988973118915344195337,
                                                     +1.29888362714502619258,
                                                     +1.21642871598462436586,
                                                     -1.22708085895067061070,
                                                     -2.03140778260244983704,
                                                     -1.6983261840450686051,
                                                     -1.6983261840450686051,
                                                     -2.03140778260244983704,
                                                     -1.22708085895067061070,
                                                     +1.21642871598462436586,
                                                     +1.29888362714502619258,
                                                     +0.988973118915344195337,
                                                     +1.43131625920334216684,
                                                     +0.521213104349852132174});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8B() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/8,
      /*time_reversible=*/true,
      /*evaluations=*/15,
      /*first_same_as_last=*/true> const integrator({+1.48819229202921310080,
                                                     -2.33864815101041943098,
                                                     +2.89105148972198900311,
                                                     -2.89688250330423987105,
                                                     +0.00378039588362668223674,
                                                     +2.89195744315817391244,
                                                     -0.00169248587771706559145,
                                                     -3.0755169612012526619,
                                                     -0.00169248587771706559145,
                                                     +2.89195744315817391244,
                                                     +0.00378039588362668223674,
                                                     -2.89688250330423987105,
                                                     +2.89105148972198900311,
                                                     -2.33864815101041943098,
                                                     +1.48819229202921310080,
                                                     +0.0},
                                                    {+0.744096146014606550401,
                                                     -0.42522792949060316509,
                                                     +0.27620166935578478606,
                                                     -0.00291550679112543397,
                                                     -1.44655105371030659441,
                                                     +1.44786891952090029734,
                                                     +1.44513247864022842343,
                                                     -1.5386047235394848638,
                                                     -1.5386047235394848638,
                                                     +1.44513247864022842343,
                                                     +1.44786891952090029734,
                                                     -1.44655105371030659441,
                                                     -0.00291550679112543397,
                                                     +0.27620166935578478606,
                                                     -0.42522792949060316509,
                                                     +0.744096146014606550401});
return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8C() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/8,
      /*time_reversible=*/true,
      /*evaluations=*/15,
      /*first_same_as_last=*/true> const integrator({+0.629030650210427818049,
                                                     +1.36934946416874222370,
                                                     -1.06458714789183904181,
                                                     +1.66335809963311356298,
                                                     -1.67896928259637402925,
                                                     -1.55946803821449795876,
                                                     +0.311790812418431890510,
                                                     +1.6589908845439910692,
                                                     +0.311790812418431890510,
                                                     -1.55946803821449795876,
                                                     -1.67896928259637402925,
                                                     +1.66335809963311356298,
                                                     -1.06458714789183904181,
                                                     +1.36934946416874222370,
                                                     +0.629030650210427818049,
                                                     +0.0},
                                                    {+0.314515325105213909024,
                                                     +0.999190057189585020872,
                                                     +0.15238115813845159094,
                                                     +0.29938547587063726059,
                                                     -0.00780559148163023314,
                                                     -1.61921866040543599400,
                                                     -0.623838612898033034124,
                                                     +0.98539084848121147984,
                                                     +0.98539084848121147984,
                                                     -0.623838612898033034124,
                                                     -1.61921866040543599400,
                                                     -0.00780559148163023314,
                                                     +0.29938547587063726059,
                                                     +0.15238115813845159094,
                                                     +0.999190057189585020872,
                                                     +0.314515325105213909024});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8D() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/8,
      /*time_reversible=*/true,
      /*evaluations=*/15,
      /*first_same_as_last=*/true> const integrator({+0.914844246229642658287,
                                                     +0.253693336566286009974,
                                                     -1.44485223686030647660,
                                                     -0.158240635368502468458,
                                                     +1.93813913762291232471,
                                                     -1.96061023297558163691,
                                                     +0.102799849392219431139,
                                                     +1.7084530707866603157,
                                                     +0.102799849392219431139,
                                                     -1.96061023297558163691,
                                                     +1.93813913762291232471,
                                                     -0.158240635368502468458,
                                                     -1.44485223686030647660,
                                                     +0.253693336566286009974,
                                                     +0.914844246229642658287,
                                                     +0.0},
                                                    {+0.457422123114821329143,
                                                     +0.584268791397964334130,
                                                     -0.595579450147010233314,
                                                     -0.801546436114404472530,
                                                     +0.88994925112720492813,
                                                     -0.01123554767633465610,
                                                     -0.92890519179168110289,
                                                     +0.90562646008943987343,
                                                     +0.90562646008943987343,
                                                     -0.92890519179168110289,
                                                     -0.01123554767633465610,
                                                     +0.88994925112720492813,
                                                     -0.801546436114404472530,
                                                     -0.595579450147010233314,
                                                     +0.584268791397964334130,
                                                     +0.457422123114821329143});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8E() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/8,
      /*time_reversible=*/true,
      /*evaluations=*/15,
      /*first_same_as_last=*/true> const integrator({+1.30300165757516838484,
                                                     +0.107990467718098279648,
                                                     -2.04809795883490205633,
                                                     +0.00536018921375238082832,
                                                     -0.0719180053650705075005,
                                                     +2.52778927318028339169,
                                                     +0.0227738840126312259937,
                                                     -2.6937990149999221983,
                                                     +0.0227738840126312259937,
                                                     +2.52778927318028339169,
                                                     -0.0719180053650705075005,
                                                     +0.00536018921375238082832,
                                                     -2.04809795883490205633,
                                                     +0.107990467718098279648,
                                                     +1.30300165757516838484,
                                                     +0.0},
                                                    {+0.651500828787584192418,
                                                     +0.705496062646633332241,
                                                     -0.97005374555840188834,
                                                     -1.02136888481057483775,
                                                     -0.0332789080756590633361,
                                                     +1.22793563390760644210,
                                                     +1.27528157859645730884,
                                                     -1.33551256549364548617,
                                                     -1.33551256549364548617,
                                                     +1.27528157859645730884,
                                                     +1.22793563390760644210,
                                                     -0.0332789080756590633361,
                                                     -1.02136888481057483775,
                                                     -0.97005374555840188834,
                                                     +0.705496062646633332241,
                                                     +0.651500828787584192418});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/3,
                                          /*first_same_as_last=*/true> const&
CandyRozmus1991ForestRuth1990() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/4,
      /*time_reversible=*/true,
      /*evaluations=*/3,
      /*first_same_as_last=*/true> const integrator({+1.3512071919596576340,
                                                     -1.7024143839193152681,
                                                     +1.3512071919596576340,
                                                     +0.0},
                                                    {+0.6756035959798288170,
                                                     -0.1756035959798288170,
                                                     -0.1756035959798288170,
                                                     +0.6756035959798288170});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/2,
                                          /*time_reversible=*/false,
                                          /*evaluations=*/2,
                                          /*first_same_as_last=*/false> const&
McLachlanAtela1992Order2Optimal() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/2,
      /*time_reversible=*/false,
      /*evaluations=*/2,
      /*first_same_as_last=*/false> const integrator({+0.7071067811865475244,
                                                      +0.2928932188134524756},
                                                     {+0.2928932188134524756,
                                                      +0.7071067811865475244});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/3,
                                          /*time_reversible=*/false,
                                          /*evaluations=*/3,
                                          /*first_same_as_last=*/false> const&
McLachlanAtela1992Order3Optimal() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/3,
      /*time_reversible=*/false,
      /*evaluations=*/3,
      /*first_same_as_last=*/false> const integrator({+0.9196615230173998571,
                                                      -0.1879916187991597820,
                                                      +0.2683300957817599250},
                                                     {+0.2683300957817599250,
                                                      -0.1879916187991597820,
                                                      +0.9196615230173998571});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/2,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/2,
                                          /*first_same_as_last=*/true> const&
McLachlan1995S2() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/2,
      /*time_reversible=*/true,
      /*evaluations=*/2,
      /*first_same_as_last=*/true> const integrator({0.5, 0.5, 0.0},
                                                    {+0.19318332750378357396,
                                                     +0.61363334499243285207,
                                                     +0.19318332750378357396});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/5,
                                          /*first_same_as_last=*/true> const&
McLachlan1995SS5() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/4,
      /*time_reversible=*/true,
      /*evaluations=*/5,
      /*first_same_as_last=*/true> const integrator({+0.28,
                                                     +0.6254664284676700450,
                                                     -0.8109328569353400900,
                                                     +0.6254664284676700450,
                                                     +0.28,
                                                     +0.0},
                                                    {+0.14,
                                                     +0.45273321423383502251,
                                                     -0.0927332142338350225,
                                                     -0.0927332142338350225,
                                                     +0.45273321423383502251,
                                                     +0.14});
    return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/4,
                                          /*first_same_as_last=*/true> const&
McLachlan1995S4() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/4,
      /*time_reversible=*/true,
      /*evaluations=*/4,
      /*first_same_as_last=*/true> const integrator({+0.54545454545454545455,
                                                     -0.045454545454545454545,
                                                     -0.045454545454545454545,
                                                     +0.54545454545454545455,
                                                     +0.0},
                                                    {+0.16913927992207204518,
                                                     -0.29918620390405079951,
                                                     +1.2600938479639575087,
                                                     -0.29918620390405079951,
                                                     +0.16913927992207204518});
  return integrator;
}
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/5,
                                          /*first_same_as_last=*/true> const&
McLachlan1995S5() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/4,
      /*time_reversible=*/true,
      /*evaluations=*/5,
      /*first_same_as_last=*/true> const
      integrator({+0.4, -0.1, 0.4, -0.1, 0.4, 0.0},
                 {+0.089269454226475244887,
                  -0.097336042636895508015,
                  +0.50806658841042026313,
                  +0.50806658841042026313,
                  -0.097336042636895508015,
                  +0.089269454226475244887});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/9,
                                          /*first_same_as_last=*/true> const&
McLachlan1995SS9() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/6,
      /*time_reversible=*/true,
      /*evaluations=*/9,
      /*first_same_as_last=*/true> const integrator({+0.1867,
                                                     +0.55549702371247839916,
                                                     +0.12946694891347535806,
                                                     -0.84326562338773460855,
                                                     +0.9432033015235617027,
                                                     -0.84326562338773460855,
                                                     +0.12946694891347535806,
                                                     +0.55549702371247839916,
                                                     +0.1867,
                                                     +0.0},
                                                    {+0.09335,
                                                     +0.37109851185623919958,
                                                     +0.34248198631297687861,
                                                     -0.35689933723712962525,
                                                     +0.0499688390679135471,
                                                     +0.0499688390679135471,
                                                     -0.35689933723712962525,
                                                     +0.34248198631297687861,
                                                     +0.37109851185623919958,
                                                     +0.09335});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
McLachlan1995SS15() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/8,
      /*time_reversible=*/true,
      /*evaluations=*/15,
      /*first_same_as_last=*/true> const integrator({+0.7416703643506129534,
                                                     -0.4091008258000315940,
                                                     +0.1907547102962383800,
                                                     -0.5738624711160822667,
                                                     +0.2990641813036559238,
                                                     +0.3346249182452981838,
                                                     +0.3152930923967665966,
                                                     -0.796887939352916354,
                                                     +0.3152930923967665966,
                                                     +0.3346249182452981838,
                                                     +0.2990641813036559238,
                                                     -0.5738624711160822667,
                                                     +0.1907547102962383800,
                                                     -0.4091008258000315940,
                                                     +0.7416703643506129534,
                                                     +0.0},
                                                    {+0.37083518217530647672,
                                                     +0.1662847692752906797,
                                                     -0.1091730577518966070,
                                                     -0.1915538804099219434,
                                                     -0.13739914490621317141,
                                                     +0.3168445497744770538,
                                                     +0.32495900532103239021,
                                                     -0.2407974234780748787,
                                                     -0.2407974234780748787,
                                                     +0.32495900532103239021,
                                                     +0.3168445497744770538,
                                                     -0.13739914490621317141,
                                                     -0.1915538804099219434,
                                                     -0.1091730577518966070,
                                                     +0.1662847692752906797,
                                                     +0.37083518217530647672});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/17,
                                          /*first_same_as_last=*/true> const&
McLachlan1995SS17() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/8,
      /*time_reversible=*/true,
      /*evaluations=*/17,
      /*first_same_as_last=*/true> const integrator({+0.12886597938144329897,
                                                     +0.5815140871052509624,
                                                     -0.41017537146985013753,
                                                     +0.1851469357165877327,
                                                     -0.40955234342085141934,
                                                     +0.14440594108001204106,
                                                     +0.27833550039367965131,
                                                     +0.31495668391629485789,
                                                     -0.626994825405133975,
                                                     +0.31495668391629485789,
                                                     +0.27833550039367965131,
                                                     +0.14440594108001204106,
                                                     -0.40955234342085141934,
                                                     +0.1851469357165877327,
                                                     -0.41017537146985013753,
                                                     +0.5815140871052509624,
                                                     +0.12886597938144329897,
                                                     +0.0},
                                                    {+0.064432989690721649485,
                                                     +0.35519003324334713070,
                                                     +0.0856693578177004124,
                                                     -0.11251421787663120244,
                                                     -0.1122027038521318433,
                                                     -0.13257320117041968914,
                                                     +0.2113707207368458462,
                                                     +0.2966460921549872546,
                                                     -0.1560190707444195585,
                                                     -0.1560190707444195585,
                                                     +0.2966460921549872546,
                                                     +0.2113707207368458462,
                                                     -0.13257320117041968914,
                                                     -0.1122027038521318433,
                                                     -0.11251421787663120244,
                                                     +0.0856693578177004124,
                                                     +0.35519003324334713070,
                                                     +0.064432989690721649485});
  return integrator;
}

template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/6,
                                          /*first_same_as_last=*/true> const&
BlanesMoan2002S6() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/4,
      /*time_reversible=*/true,
      /*evaluations=*/6,
      /*first_same_as_last=*/true> const integrator({+0.20951510661336200,
                                                     -0.14385177317981800,
                                                     +0.43433666656645600,
                                                     +0.43433666656645600,
                                                     -0.14385177317981800,
                                                     +0.20951510661336200,
                                                     +0.0},
                                                    {+0.079203696431195700,
                                                     +0.35317290604977400,
                                                     -0.042065080357719500,
                                                     +0.21937695575349960,
                                                     -0.042065080357719500,
                                                     +0.35317290604977400,
                                                     +0.079203696431195700});
  return integrator;
}
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/10,
                                          /*first_same_as_last=*/true> const&
BlanesMoan2002S10() {
  static SymplecticPartitionedRungeKuttaIntegrator<
      Position,
      Momentum,
      /*order=*/6,
      /*time_reversible=*/true,
      /*evaluations=*/10,
      /*first_same_as_last=*/true> const integrator({+0.14881644790104200,
                                                     -0.13238586576778400,
                                                     +0.067307604692185000,
                                                     +0.43266640257817500,
                                                     -0.01640458940361800,
                                                     -0.01640458940361800,
                                                     +0.43266640257817500,
                                                     +0.067307604692185000,
                                                     -0.13238586576778400,
                                                     +0.14881644790104200,
                                                     +0.0},
                                                    {+0.050262764400392200,
                                                     +0.41351430042834400,
                                                     +0.045079889794397700,
                                                     -0.18805485381956900,
                                                     +0.54196067845078000,
                                                     -0.7255255585086898,
                                                     +0.54196067845078000,
                                                     -0.18805485381956900,
                                                     +0.045079889794397700,
                                                     +0.41351430042834400,
                                                     +0.050262764400392200});
  return integrator;
}

}  // namespace integrators
}  // namespace principia
