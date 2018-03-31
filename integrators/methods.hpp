
#pragma once

#include "base/mod.hpp"
#include "base/not_constructible.hpp"
#include "numerics/fixed_arrays.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace methods {

using base::mod;
using base::not_constructible;
using numerics::FixedStrictlyLowerTriangularMatrix;
using numerics::FixedVector;

struct EmbeddedExplicitRungeKuttaNyström : not_constructible {
  // static constexpr int higher_order = ...;
  // static constexpr int lower_order = ...;
  // static constexpr int stages = ...;
  // static constexpr bool first_same_as_last = ...;
  // static constexpr serialization::AdaptiveStepSizeIntegrator::Kind kind = ..;
  // static constexpr FixedVector<double, stages> c = ...;
  // static constexpr FixedStrictlyLowerTriangularMatrix<double, stages> a = ..;
  // static constexpr FixedVector<double, stages> b_hat = ...;
  // static constexpr FixedVector<double, stages> b_prime_hat = ...;
  // static constexpr FixedVector<double, stages> b = ...;
  // static constexpr FixedVector<double, stages> b_prime = ...;
};

struct SymmetricLinearMultistep : not_constructible {
  static constexpr int Half(int const order) {
    return order / 2 + 1;
  }
  // static constexpr int order = ...;
  // static constexpr serialization::FixedStepSizeIntegrator::Kind kind = ...;
  // static constexpr FixedVector<double, Half(order)> const ɑ(...);
  // static constexpr FixedVector<double, Half(order)> const β_numerator(...);
  // static constexpr double β_denominator = ...;
};

struct SymplecticRungeKuttaNyström : not_constructible {
 protected:
  using CompositionMethod =
      serialization::FixedStepSizeIntegrator::CompositionMethod;
  static constexpr CompositionMethod BA =
      serialization::FixedStepSizeIntegrator::BA;
  static constexpr CompositionMethod ABA =
      serialization::FixedStepSizeIntegrator::ABA;
  static constexpr CompositionMethod BAB =
      serialization::FixedStepSizeIntegrator::BAB;

  static constexpr int Stages(int const evaluations,
                              CompositionMethod const composition) {
    return composition == BA ? evaluations : evaluations + 1;
  }
  // static constexpr int order = ...;
  // static constexpr bool time_reversible = ...;
  // static constexpr int evaluations = ...;
  // static constexpr CompositionMethod composition = ...;
  // static constexpr serialization::FixedStepSizeIntegrator::Kind kind = ...;
  // static constexpr int stages = Stages(evaluations, composition);
  // static constexpr FixedVector<double, stages> a(...);
  // static constexpr FixedVector<double, stages> b(...);
};

struct SymplecticPartitionedRungeKutta : not_constructible {
  static constexpr int Stages(int const evaluations,
                              bool const first_same_as_last) {
    return first_same_as_last ? evaluations + 1 : evaluations;
  }
  // static constexpr int order = ...;
  // static constexpr bool time_reversible = ...;
  // static constexpr int evaluations = ...;
  // static constexpr bool first_same_as_last = ...;
  // static constexpr serialization::FixedStepSizeIntegrator::Kind kind = ..;
  // static constexpr int stages = Stages(evaluations, first_same_as_last);
  // static constexpr FixedVector<double, stages> a(...);
  // static constexpr FixedVector<double, stages> b(...);
};

// Every SPRK may be transformed into an SRKN by specifying a composition method
// (the possible composition methods are constrained by the properties of the
// SPRK).  This struct effects that transformation.
template<typename SymplecticPartitionedRungeKuttaMethod,
         serialization::FixedStepSizeIntegrator::CompositionMethod composition_>
struct AsSymplecticRungeKuttaNyström {
  static_assert(std::is_base_of<SymplecticPartitionedRungeKutta,
                                SymplecticPartitionedRungeKuttaMethod>::value,
                "Method must be derived from SymplecticPartitionedRungeKutta");
  static_assert(
      !SymplecticPartitionedRungeKuttaMethod::first_same_as_last ||
          composition_ == serialization::FixedStepSizeIntegrator::ABA ||
          composition_ == serialization::FixedStepSizeIntegrator::BAB,
      "requested |composition| must be ABA or BAB for this first-same-as-last "
      "method");
  static_assert(SymplecticPartitionedRungeKuttaMethod::first_same_as_last ||
                    composition_ == serialization::FixedStepSizeIntegrator::BA,
                "requested |composition| must be BA for this method which is "
                "not first-same-as-last");
  static_assert(composition_ != serialization::FixedStepSizeIntegrator::ABA,
                "ABA not supported until C++17");

  struct Method : SymplecticRungeKuttaNyström {
    static constexpr int order = SymplecticPartitionedRungeKuttaMethod::order;
    static constexpr bool time_reversible =
        SymplecticPartitionedRungeKuttaMethod::time_reversible;
    static constexpr int evaluations =
        SymplecticPartitionedRungeKuttaMethod::evaluations;
    static constexpr CompositionMethod composition = composition_;
    static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
        SymplecticPartitionedRungeKuttaMethod::kind;
    static constexpr int stages = Stages(evaluations, composition);

    static constexpr FixedVector<double, stages> Shift(
        FixedVector<double, stages> const& a,
        CompositionMethod const composition) {
#if 0
      if (composition == ABA) {
        FixedVector<double, stages> shifted_a;
        // |*this| is a |BAB| method, with A and B interchangeable.  Exchanging
        // A and B shifts |a_| (because |ABA| means b₀ vanishes, whereas |BAB|
        // means aᵣ vanishes).
        for (int i = 0; i < stages; ++i) {
          shifted_a[i] = a[mod(i - 1, stages)];
        }
        return shifted_a;
      } else {
#endif
        return a;
#if 0
      }
#endif
  }

    static constexpr FixedVector<double, stages> a =
        Shift(SymplecticPartitionedRungeKuttaMethod::a, composition);
    static constexpr FixedVector<double, stages> b =
        SymplecticPartitionedRungeKuttaMethod::b;
  };
};


// The following methods have coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
struct BlanesMoan2002S6 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 4;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 6;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::BLANES_MOAN_2002_S6;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.20951510661336200,
                                                   -0.14385177317981800,
                                                   +0.43433666656645600,
                                                   +0.43433666656645600,
                                                   -0.14385177317981800,
                                                   +0.20951510661336200,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.079203696431195700,
                                                   +0.35317290604977400,
                                                   -0.042065080357719500,
                                                   +0.21937695575349960,
                                                   -0.042065080357719500,
                                                   +0.35317290604977400,
                                                   +0.079203696431195700}}};
};
struct BlanesMoan2002S10 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 6;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 10;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::BLANES_MOAN_2002_S10;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.14881644790104200,
                                                   -0.13238586576778400,
                                                   +0.067307604692185000,
                                                   +0.43266640257817500,
                                                   -0.01640458940361800,
                                                   -0.01640458940361800,
                                                   +0.43266640257817500,
                                                   +0.067307604692185000,
                                                   -0.13238586576778400,
                                                   +0.14881644790104200,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.050262764400392200,
                                                   +0.41351430042834400,
                                                   +0.045079889794397700,
                                                   -0.18805485381956900,
                                                   +0.54196067845078000,
                                                   -0.7255255585086898,
                                                   +0.54196067845078000,
                                                   -0.18805485381956900,
                                                   +0.045079889794397700,
                                                   +0.41351430042834400,
                                                   +0.050262764400392200}}};
};
struct BlanesMoan2002SRKN6B : SymplecticRungeKuttaNyström {
  static constexpr int order = 4;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 6;
  static constexpr CompositionMethod composition = BAB;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::BLANES_MOAN_2002_SRKN_6B;
  static constexpr int stages = Stages(evaluations, composition);
  static constexpr FixedVector<double, stages> a{{{0.24529895718427100,
                                                   0.60487266571108000,
                                                   -0.35017162289535100,
                                                   -0.35017162289535100,
                                                   0.60487266571108000,
                                                   0.24529895718427100,
                                                   0.0}}};
  static constexpr FixedVector<double, stages> b{{{0.082984406417405200,
                                                   0.39630980149836800,
                                                   -0.039056304922348600,
                                                   0.1195241940131508,
                                                   -0.039056304922348600,
                                                   0.39630980149836800,
                                                   0.082984406417405200}}};
};
struct BlanesMoan2002SRKN11B : SymplecticRungeKuttaNyström {
  static constexpr int order = 6;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 11;
  static constexpr CompositionMethod composition = BAB;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::BLANES_MOAN_2002_SRKN_11B;
  static constexpr int stages = Stages(evaluations, composition);
  static constexpr FixedVector<double, stages> a{{{0.12322977594627100,
                                                   0.29055379779955800,
                                                   -0.12704921262541700,
                                                   -0.24633176106207500,
                                                   0.35720887279592800,
                                                   0.2047770542914700,
                                                   0.35720887279592800,
                                                   -0.24633176106207500,
                                                   -0.12704921262541700,
                                                   0.29055379779955800,
                                                   0.12322977594627100,
                                                   0.0}}};
  static constexpr FixedVector<double, stages> b{{{0.041464998518262400,
                                                   0.19812867191806700,
                                                   -0.040006192104153300,
                                                   0.075253984301580700,
                                                   -0.011511387420687900,
                                                   0.23666992478693110,
                                                   0.23666992478693110,
                                                   -0.011511387420687900,
                                                   0.075253984301580700,
                                                   -0.040006192104153300,
                                                   0.19812867191806700,
                                                   0.041464998518262400}}};
};
struct BlanesMoan2002SRKN14A : SymplecticRungeKuttaNyström {
  static constexpr int order = 6;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 14;
  static constexpr CompositionMethod composition = ABA;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::BLANES_MOAN_2002_SRKN_14A;
  static constexpr int stages = Stages(evaluations, composition);
  static constexpr FixedVector<double, stages> a{{{0.037859319840611600,
                                                   0.10263563310243500,
                                                   -0.025867888266558700,
                                                   0.31424140307144700,
                                                   -0.13014445951741500,
                                                   0.10641770036954300,
                                                   -0.0087942431285105800,
                                                   0.2073050690568954,
                                                   -0.0087942431285105800,
                                                   0.10641770036954300,
                                                   -0.13014445951741500,
                                                   0.31424140307144700,
                                                   -0.025867888266558700,
                                                   0.10263563310243500,
                                                   0.037859319840611600}}};
  static constexpr FixedVector<double, stages> b{{{0.0,
                                                   0.091719152624461650,
                                                   0.18398317000500600,
                                                   -0.056534365832888270,
                                                   0.0049146887747128540,
                                                   0.14376112716835800,
                                                   0.32856769374680400,
                                                   -0.19641146648645423,
                                                   -0.19641146648645423,
                                                   0.32856769374680400,
                                                   0.14376112716835800,
                                                   0.0049146887747128540,
                                                   -0.056534365832888270,
                                                   0.18398317000500600,
                                                   0.091719152624461650}}};
};

// Coefficients from Forest and Ruth (1990),
// Fourth-order symplectic integration, equation 4.8.
// http://zwe.web.cern.ch/zwe/CAS/biblio/ruth-forest.pdf.
// This scheme was independently discovered by Candy and Rozmus (1991),
// A Symplectic Integration Algorithm for Separable Hamiltonian Functions
// (submitted earlier and published later than the Forest and Ruth paper).
struct CandyRozmus1991ForestRuth1990 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 4;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 3;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::
          CANDY_ROZMUS_1991_FOREST_RUTH_1990;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+1.3512071919596576340,
                                                   -1.7024143839193152681,
                                                   +1.3512071919596576340,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.6756035959798288170,
                                                   -0.1756035959798288170,
                                                   -0.1756035959798288170,
                                                   +0.6756035959798288170}}};
};

// Coefficients from Dormand, El-Mikkawy and Prince (1986),
// Families of Runge-Kutta-Nyström formulae, table 3 (the RK4(3)4FM).
// Minimizes the 4th order truncation error.
struct DormandElMikkawyPrince1986RKN434FM :
    EmbeddedExplicitRungeKuttaNyström {
  static constexpr int higher_order = 4;
  static constexpr int lower_order = 3;
  static constexpr int stages = 4;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::AdaptiveStepSizeIntegrator::Kind kind =
      serialization::AdaptiveStepSizeIntegrator::
          DORMAND_ELMIKKAWY_PRINCE_1986_RKN_434FM;
  static constexpr FixedVector<double, stages> c{{
      { 0.0         ,   1.0 /   4.0,   7.0 /  10.0,  1.0}}};
  static constexpr FixedStrictlyLowerTriangularMatrix<double, stages> a{{
      { 1.0 /   32.0,
        7.0 / 1000.0, 119.0 / 500.0,
        1.0 /   14.0,   8.0 /  27.0,  25.0 / 189.0}}};
  static constexpr FixedVector<double, stages> b_hat{{
      { 1.0 /   14.0,   8.0 /  27.0,  25.0 / 189.0,  0.0}}};
  static constexpr FixedVector<double, stages> b_prime_hat{{
      { 1.0 /   14.0,  32.0 /  81.0, 250.0 / 567.0,  5.0 / 54.0}}};
  static constexpr FixedVector<double, stages> b{{
      {-7.0 /  150.0,  67.0 / 150.0,   3.0 /  20.0, -1.0 / 20.0}}};
  static constexpr FixedVector<double, stages> b_prime{{
      {13.0 /   21.0, -20.0 /  27.0, 275.0 / 189.0, -1.0 /  3.0}}};
};

// The following methods have coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
struct McLachlan1995S2 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 2;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 2;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_S2;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{0.5, 0.5, 0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.19318332750378357396,
                                                   +0.61363334499243285207,
                                                   +0.19318332750378357396}}};
};
struct McLachlan1995S4 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 4;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 4;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_S4;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.54545454545454545455,
                                                   -0.045454545454545454545,
                                                   -0.045454545454545454545,
                                                   +0.54545454545454545455,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.16913927992207204518,
                                                   -0.29918620390405079951,
                                                   +1.2600938479639575087,
                                                   -0.29918620390405079951,
                                                   +0.16913927992207204518}}};
};
struct McLachlan1995S5 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 4;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 5;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_S5;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{
      {{+0.4, -0.1, 0.4, -0.1, 0.4, 0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.089269454226475244887,
                                                   -0.097336042636895508015,
                                                   +0.50806658841042026313,
                                                   +0.50806658841042026313,
                                                   -0.097336042636895508015,
                                                   +0.089269454226475244887}}};
};
struct McLachlan1995SB3A4 : SymplecticRungeKuttaNyström {
  static constexpr int order = 4;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 4;
  static constexpr CompositionMethod composition = ABA;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_SB3A_4;
  static constexpr int stages = Stages(evaluations, composition);
  static constexpr FixedVector<double, stages> a{{{0.18819521776883821787,
                                                   -0.021528551102171551201,
                                                   0.66666666666666666667,
                                                   -0.021528551102171551201,
                                                   0.18819521776883821787}}};
  static constexpr FixedVector<double, stages> b{{{0.0, 1.0, -0.5, -0.5, 1.0}}};
};
struct McLachlan1995SB3A5 : SymplecticRungeKuttaNyström {
  static constexpr int order = 4;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 5;
  static constexpr CompositionMethod composition = ABA;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_SB3A_5;
  static constexpr int stages = Stages(evaluations, composition);
  static constexpr FixedVector<double, stages> a{{{0.4051886183952522772,
                                                   -0.2871440408165240890,
                                                   0.3819554224212718118,
                                                   0.3819554224212718118,
                                                   -0.2871440408165240890,
                                                   0.4051886183952522772}}};
  static constexpr FixedVector<double, stages> b{{{0.0,
                                                   -0.041095890410958904110,
                                                   0.28813559322033898305,
                                                   0.50592059438123984212,
                                                   0.28813559322033898305,
                                                   -0.041095890410958904110}}};
};
struct McLachlan1995SS5 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 4;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 5;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_SS5;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.28,
                                                   +0.6254664284676700450,
                                                   -0.8109328569353400900,
                                                   +0.6254664284676700450,
                                                   +0.28,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.14,
                                                   +0.45273321423383502251,
                                                   -0.0927332142338350225,
                                                   -0.0927332142338350225,
                                                   +0.45273321423383502251,
                                                   +0.14}}};
};
struct McLachlan1995SS9 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 6;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 9;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_SS9;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.1867,
                                                   +0.55549702371247839916,
                                                   +0.12946694891347535806,
                                                   -0.84326562338773460855,
                                                   +0.9432033015235617027,
                                                   -0.84326562338773460855,
                                                   +0.12946694891347535806,
                                                   +0.55549702371247839916,
                                                   +0.1867,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.09335,
                                                   +0.37109851185623919958,
                                                   +0.34248198631297687861,
                                                   -0.35689933723712962525,
                                                   +0.0499688390679135471,
                                                   +0.0499688390679135471,
                                                   -0.35689933723712962525,
                                                   +0.34248198631297687861,
                                                   +0.37109851185623919958,
                                                   +0.09335}}};
};
struct McLachlan1995SS15 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 8;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 15;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_SS15;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.7416703643506129534,
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
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.37083518217530647672,
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
                                                   +0.37083518217530647672}}};
};
struct McLachlan1995SS17 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 8;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 17;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_SS17;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.12886597938144329897,
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
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.064432989690721649485,
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
                                                   +0.064432989690721649485}}};
};

// The following methods have coefficients from Robert I. McLachlan and Pau
// Atela (1992), The accuracy of symplectic integrators, table 2.
// http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
struct McLachlanAtela1992Order2Optimal : SymplecticPartitionedRungeKutta {
  static constexpr int order = 2;
  static constexpr bool time_reversible = false;
  static constexpr int evaluations = 2;
  static constexpr bool first_same_as_last = false;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::
          MCLACHLAN_ATELA_1992_ORDER_2_OPTIMAL;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{
      {{+0.7071067811865475244, +0.2928932188134524756}}};
  static constexpr FixedVector<double, stages> b{
      {{+0.2928932188134524756, +0.7071067811865475244}}};
};
struct McLachlanAtela1992Order3Optimal : SymplecticPartitionedRungeKutta {
  static constexpr int order = 3;
  static constexpr bool time_reversible = false;
  static constexpr int evaluations = 3;
  static constexpr bool first_same_as_last = false;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::
          MCLACHLAN_ATELA_1992_ORDER_3_OPTIMAL;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.9196615230173998571,
                                                   -0.1879916187991597820,
                                                   +0.2683300957817599250}}};
  static constexpr FixedVector<double, stages> b{{{+0.2683300957817599250,
                                                   -0.1879916187991597820,
                                                   +0.9196615230173998571}}};
};
struct McLachlanAtela1992Order4Optimal : SymplecticRungeKuttaNyström {
  static constexpr int order = 4;
  static constexpr bool time_reversible = false;
  static constexpr int evaluations = 4;
  static constexpr CompositionMethod composition = BA;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::
          MCLACHLAN_ATELA_1992_ORDER_4_OPTIMAL;
  static constexpr int stages = Stages(evaluations, composition);
  static constexpr FixedVector<double, stages> a{{{0.5153528374311229364,
                                                   -0.085782019412973646,
                                                   0.4415830236164665242,
                                                   0.1288461583653841854}}};
  static constexpr FixedVector<double, stages> b{{{0.1344961992774310892,
                                                   -0.2248198030794208058,
                                                   0.7563200005156682911,
                                                   0.3340036032863214255}}};
};
struct McLachlanAtela1992Order5Optimal : SymplecticRungeKuttaNyström {
  static constexpr int order = 5;
  static constexpr bool time_reversible = false;
  static constexpr int evaluations = 6;
  static constexpr CompositionMethod composition = BA;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::
          MCLACHLAN_ATELA_1992_ORDER_5_OPTIMAL;
  static constexpr int stages = Stages(evaluations, composition);
  static constexpr FixedVector<double, stages> a{{{0.339839625839110000,
                                                   -0.088601336903027329,
                                                   0.5858564768259621188,
                                                   -0.603039356536491888,
                                                   0.3235807965546976394,
                                                   0.4423637942197494587}}};
  static constexpr FixedVector<double, stages> b{{{0.1193900292875672758,
                                                   0.6989273703824752308,
                                                   -0.1713123582716007754,
                                                   0.4012695022513534480,
                                                   0.0107050818482359840,
                                                   -0.0589796254980311632}}};
};

// This integrator goes by many names, see Hairer, Lubich, and Wanner (2003),
// Geometric numerical integration illustrated by the Störmer–Verlet method,
// http://www.math.kit.edu/ianm3/lehre/geonumint2009s/media/gni_by_stoermer-verlet.pdf,
// section 1.3 (Historical remarks).
// Notably, it appears in Philosophiae Naturalis Principia Mathematica,
// in section II (De inventione Virium centripetarum), in the proof of
// theorem I.  See p. 37 of the first edition,
// https://cudl.lib.cam.ac.uk/view/PR-ADV-B-00039-00001/97.
// It also appears in:
// - Delambre (1790), De l'usage du calcul differentiel dans la
//   construction des tables astronomiques, in Mémoires de l'Académie Royale des
//   Sciences de Turin, vol. V, Mémoires présentés à l'Académie, p. 143-180.
//   http://www.biodiversitylibrary.org/item/32318#page/698/mode/1up.
// - Størmer (1907), Sur les trajectoires des corpuscules électrisés dans
//   l'espace, avec application aux aurores boréales.
//   https://hal.archives-ouvertes.fr/jpa-00242574/document.
// - Verlet (1967) Computer "Experiments" on classical fluids. I.
//   Thermodynamical properties of Lennard-Jones molecules.
//   http://www.chemie.unibas.ch/~steinhauser/teaching/FS2014/MD-Simulation/Handout_12.pdf.
struct NewtonDelambreStørmerVerletLeapfrog : SymplecticPartitionedRungeKutta {
  static constexpr int order = 2;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 1;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::
          NEWTON_DELAMBRE_STORMER_VERLET_LEAPFROG;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{1.0, 0.0}}};
  static constexpr FixedVector<double, stages> b{{{0.5, 0.5}}};
};

// Coefficients from Okunbor and Skeel (1994),
// Canonical Runge-Kutta-Nyström methods of orders 5 and 6,
// http://bionum.cs.purdue.edu/94OkSk.pdf.
// NOTE(egg): The coefficients were actually copied from McLachlan (1995), they
// seem to differ after a dozen significant figures or so.  Okunbor and Skeel
// remark "we did not use HYBRJ1 to improve the accuracy of method coefficients
// as we did in section 3.1".  We assume McLachlan's version is accurate.
// TODO(egg): Derive the coefficients with Mathematica.
struct OkunborSkeel1994Order6Method13 : SymplecticRungeKuttaNyström {
  static constexpr int order = 6;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 7;
  static constexpr CompositionMethod composition = ABA;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::
          OKUNBOR_SKEEL_1994_ORDER_6_METHOD_13;
  static constexpr int stages = Stages(evaluations, composition);
  static constexpr FixedVector<double, stages> a{{{-1.0130879789171747298,
                                                   1.1874295737325427070,
                                                   -0.018335852096460590340,
                                                   0.34399425728109261313,
                                                   0.34399425728109261313,
                                                   -0.018335852096460590340,
                                                   1.1874295737325427070,
                                                   -1.0130879789171747298}}};
  static constexpr FixedVector<double, stages> b{{{0.0,
                                                   0.00016600692650009894,
                                                   -0.37962421426377360608,
                                                   0.68913741185181063674,
                                                   0.38064159097092574080,
                                                   0.68913741185181063674,
                                                   -0.37962421426377360608,
                                                   0.00016600692650009894}}};
};

// The following methods are from Quinlan (1999),
// Resonances and instabilities in symmetric multistep methods,
// https://arxiv.org/abs/astro-ph/9901136.
struct Quinlan1999Order8A : SymmetricLinearMultistep {
  static constexpr int order = 8;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8A;
  static constexpr FixedVector<double, Half(order)> const ɑ{
      {{1.0, -2.0, 2.0, -2.0, 2.0}}};
  static constexpr FixedVector<double, Half(order)> const β_numerator{
      {{0.0, 22081.0, -29418.0, 75183.0, -75212.0}}};
  static constexpr double β_denominator = 15120.0;
};
struct Quinlan1999Order8B : SymmetricLinearMultistep {
  static constexpr int order = 8;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8B;
  static constexpr FixedVector<double, Half(order)> const ɑ{
      {{1.0, 0.0, 0.0, -0.5, -1.0}}};
  static constexpr FixedVector<double, Half(order)> const β_numerator{
      {{0.0, 192481.0, 6582.0, 816783.0, -156812.0}}};
  static constexpr double β_denominator = 120960.0;
};

// The following methods are from Quinlan and Tremaine (1990),
// Symmetric multistep methods for the numerical integration of planetary
// orbits,
// http://adsabs.harvard.edu/full/1990AJ....100.1694Q.
struct QuinlanTremaine1990Order8 : SymmetricLinearMultistep {
  static constexpr int order = 8;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_8;
  static constexpr FixedVector<double, Half(order)> const ɑ{
      {{1.0, -2.0, 2.0, -1.0, 0.0}}};
  static constexpr FixedVector<double, Half(order)> const β_numerator{
      {{0.0, 17671.0, -23622.0, 61449.0, -50516.0}}};
  static constexpr double β_denominator = 12096.0;
};
struct QuinlanTremaine1990Order10 : SymmetricLinearMultistep {
  static constexpr int order = 10;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_10;
  static constexpr FixedVector<double, Half(order)> const ɑ{
      {{1.0, -1.0, 1.0, -1.0, 1.0, -2.0}}};
  static constexpr FixedVector<double, Half(order)> const β_numerator{
      {{0.0, 399187.0, -485156.0, 2391436.0, -2816732.0, 4651330.0}}};
  static constexpr double β_denominator = 241920.0;
};
struct QuinlanTremaine1990Order12 : SymmetricLinearMultistep {
  static constexpr int order = 12;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_12;
  static constexpr FixedVector<double, Half(order)> const ɑ{
      {{1.0, -2.0, 2.0, -1.0, 0.0, 0.0, 0.0}}};
  static constexpr FixedVector<double, Half(order)> const β_numerator{
      {{0.0,
        90987349.0,
        -229596838.0,
        812627169.0,
        -1628539944.0,
        2714971338.0,
        -3041896548.0}}};
  static constexpr double β_denominator = 53222400.0;
};
struct QuinlanTremaine1990Order14 : SymmetricLinearMultistep {
  static constexpr int order = 14;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_14;
  static constexpr FixedVector<double, Half(order)> const ɑ{
      {{1.0, -2.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0}}};
  static constexpr FixedVector<double, Half(order)> const β_numerator{
      {{0.0,
        433489274083.0,
        -1364031998256.0,
        5583113380398.0,
        -14154444148720.0,
        28630585332045.0,
        -42056933842656.0,
        48471792742212.0}}};
  static constexpr double β_denominator = 237758976000.0;
};

// Coefficients from Ruth (1983), A canonical integration technique,
// https://accelconf.web.cern.ch/accelconf/p83/PDF/PAC1983_2669.PDF.
struct Ruth1983 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 3;
  static constexpr bool time_reversible = false;
  static constexpr int evaluations = 3;
  static constexpr bool first_same_as_last = false;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::RUTH_1983;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{
      {{2.0 / 3.0, -2.0 / 3.0, 1.0}}};
  static constexpr FixedVector<double, stages> b{
      {{7.0 / 24.0, 3.0 / 4.0, -1.0 / 24.0}}};
};

// Coefficients from Suzuki (1990), Fractal decomposition of exponential
// operators with applications to many-body theories and Monte Carlo
// simulations.
struct Suzuki1990 : SymplecticPartitionedRungeKutta {
  static constexpr int order = 4;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 5;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::SUZUKI_1990;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.41449077179437573714,
                                                   +0.41449077179437573714,
                                                   -0.65796308717750294857,
                                                   +0.41449077179437573714,
                                                   +0.41449077179437573714,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.20724538589718786857,
                                                   +0.41449077179437573714,
                                                   -0.12173615769156360571,
                                                   -0.12173615769156360571,
                                                   +0.41449077179437573714,
                                                   +0.20724538589718786857}}};
};

// The following methods have coefficients from Yoshida (1990),
// Construction of higher order symplectic integrators
// http://sixtrack.web.cern.ch/SixTrack/doc/yoshida00.pdf.
// NOTE(egg): The coefficients were derived from equations 5.4 through 5.17
// rather than computed from the wᵢ given in tables 1 and 2.  The results were
// then cross-checked against those obtained from the tables.
struct Yoshida1990Order6A : SymplecticPartitionedRungeKutta {
  static constexpr int order = 6;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 7;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::YOSHIDA_1990_ORDER_6A;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.78451361047755726382,
                                                   +0.23557321335935813369,
                                                   -1.17767998417887100695,
                                                   +1.31518632068391121889,
                                                   -1.17767998417887100695,
                                                   +0.23557321335935813369,
                                                   +0.78451361047755726382,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.392256805238778631910,
                                                   +0.51004341191845769875,
                                                   -0.47105338540975643663,
                                                   +0.06875316825252010597,
                                                   +0.06875316825252010597,
                                                   -0.47105338540975643663,
                                                   +0.51004341191845769875,
                                                   +0.392256805238778631910}}};
};
struct Yoshida1990Order6B : SymplecticPartitionedRungeKutta {
  static constexpr int order = 6;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 7;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::YOSHIDA_1990_ORDER_6B;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+1.43984816797678309102,
                                                   +0.00426068187079201616,
                                                   -2.13228522200145152088,
                                                   +2.37635274430775282740,
                                                   -2.13228522200145152088,
                                                   +0.00426068187079201616,
                                                   +1.43984816797678309102,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.71992408398839154551,
                                                   +0.72205442492378755359,
                                                   -1.06401227006532975236,
                                                   +0.12203376115315065326,
                                                   +0.12203376115315065326,
                                                   -1.06401227006532975236,
                                                   +0.72205442492378755359,
                                                   +0.71992408398839154551}}};
};
struct Yoshida1990Order6C : SymplecticPartitionedRungeKutta {
  static constexpr int order = 6;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 7;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::YOSHIDA_1990_ORDER_6C;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+1.44778256239929793290,
                                                   -2.14403531630538931060,
                                                   +0.00152886228424927025338,
                                                   +2.38944778324368421490,
                                                   +0.00152886228424927025338,
                                                   -2.14403531630538931060,
                                                   +1.44778256239929793290,
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.72389128119964896645,
                                                   -0.34812637695304568885,
                                                   -1.07125322701057002017,
                                                   +1.19548832276396674257,
                                                   +1.19548832276396674257,
                                                   -1.07125322701057002017,
                                                   -0.34812637695304568885,
                                                   +0.72389128119964896645}}};
};
struct Yoshida1990Order8A : SymplecticPartitionedRungeKutta {
  static constexpr int order = 8;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 15;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::YOSHIDA_1990_ORDER_8A;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+1.04242620869970426435,
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
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.521213104349852132174,
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
                                                   +0.521213104349852132174}}};
};
struct Yoshida1990Order8B : SymplecticPartitionedRungeKutta {
  static constexpr int order = 8;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 15;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::YOSHIDA_1990_ORDER_8B;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+1.48819229202921310080,
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
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.744096146014606550401,
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
                                                   +0.744096146014606550401}}};
};
struct Yoshida1990Order8C : SymplecticPartitionedRungeKutta {
  static constexpr int order = 8;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 15;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::YOSHIDA_1990_ORDER_8C;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.629030650210427818049,
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
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.314515325105213909024,
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
                                                   +0.314515325105213909024}}};
};
struct Yoshida1990Order8D : SymplecticPartitionedRungeKutta {
  static constexpr int order = 8;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 15;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::YOSHIDA_1990_ORDER_8D;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+0.914844246229642658287,
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
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.457422123114821329143,
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
                                                   +0.457422123114821329143}}};
};
struct Yoshida1990Order8E : SymplecticPartitionedRungeKutta {
  static constexpr int order = 8;
  static constexpr bool time_reversible = true;
  static constexpr int evaluations = 15;
  static constexpr bool first_same_as_last = true;
  static constexpr serialization::FixedStepSizeIntegrator::Kind kind =
      serialization::FixedStepSizeIntegrator::YOSHIDA_1990_ORDER_8E;
  static constexpr int stages = Stages(evaluations, first_same_as_last);
  static constexpr FixedVector<double, stages> a{{{+1.30300165757516838484,
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
                                                   +0.0}}};
  static constexpr FixedVector<double, stages> b{{{+0.651500828787584192418,
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
                                                   +0.651500828787584192418}}};
};

}  // namespace methods
}  // namespace integrators
}  // namespace principia
