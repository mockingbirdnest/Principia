
#pragma once

#include <type_traits>

#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_runge_kutta_nyström_integrator {

using numerics::FixedVector;

// A symplectic partitioned Runge-Kutta integrator.  Does not subclass
// |Integrator|; used to generate (less general)
// |SymplecticRungeKuttaNyströmIntegrator|s.
// Represents a single-step method for the solution of
//   (q, p)′ = X(q, p, t), with X = A(q, p, t) + B(q, p, t).
// |Position| is the type of |q|, and |Momentum| is that of |p|.
// The step is the composition of evolutions
//   exp(aᵣ₋₁ h A) exp(bᵣ₋₁ h B) ... exp(a₀ h A) exp(b₀ h B);
// A and B are interchangeable.  If aᵣ₋₁ vanishes, this becomes
//   exp(bᵣ₋₁ h B) exp(aᵣ₋₂ h A) ... exp(a₀ h A) exp(b₀ h B).
// In that case, |first_same_as_last| is true.
// The equation solved by this integrator is a general case of that solved by
// a |SymplecticRungeKuttaNyströmIntegrator|, see equation (1) in the
// appropriate file.  This may therefore be turned into a
// |SymplecticRungeKuttaNyströmIntegrator|.
// In the |first_same_as_last| case, since A and B are interchangeable
// for a |SymplecticPartitionedRungeKuttaIntegrator|, the step
//   exp(bᵣ₋₁ h A) exp(aᵣ₋₂ h B) ... exp(a₀ h B) exp(b₀ h A)
// also provides an integrator.  Since for a
// |SymplecticRungeKuttaNyströmIntegrator| A and B are *not* interchangeable,
// because of the requirement that [B, [B, [B, A]]] = 0 in (1), there are two
// different |SymplecticRungeKuttaNyströmIntegrator| corresponding to a
// |first_same_as_last| |SymplecticPartitionedRungeKuttaIntegrator|: one is
// |ABA|, the other is |BAB|.
// NOTE(egg): The |SymplecticRungeKuttaNyströmIntegrator| thus constructed will
// serialize as a |DUMMY| and probably break in all sorts of hilarious ways if
// deserialized.
// TODO(egg): Make them serializable/deserializable.  We need to prevent
// combinatorial explosion.
template<typename Position,
         typename Momentum,
         int order_,
         bool time_reversible_,
         int evaluations_,
         bool first_same_as_last_>
class SymplecticPartitionedRungeKuttaIntegrator {
  static constexpr int stages_ = first_same_as_last_ ? evaluations_ + 1
                                                     : evaluations_;

 public:
  SymplecticPartitionedRungeKuttaIntegrator(
      FixedVector<double, stages_> const& a,
      FixedVector<double, stages_> const& b);

  static constexpr int order = order_;
  static constexpr int evaluations = evaluations_;
  static constexpr bool time_reversible = time_reversible_;
  static constexpr bool first_same_as_last = first_same_as_last_;

  // If |first_same_as_last|, |composition_method| must be |ABA| or |BAB|.
  // Otherwise, it must be |BA|.
  template<CompositionMethod composition_method>
  SymplecticRungeKuttaNyströmIntegrator<
      Position,
      order_,
      time_reversible_,
      evaluations_,
      composition_method> const& AsRungeKuttaNyströmIntegrator() const;

 private:
  static constexpr auto BA = methods::SymplecticRungeKuttaNyström::BA;
  static constexpr auto ABA = methods::SymplecticRungeKuttaNyström::ABA;
  static constexpr auto BAB = methods::SymplecticRungeKuttaNyström::BAB;

  FixedVector<double, stages_> const a_;
  FixedVector<double, stages_> const b_;

  // The Runge-Kutta-Nyström methods are stored here, so that we can use them
  // by const-reference as we do for the others.  Since |*this| should be a
  // static object, this similarly obviates questions of lifetime.
  // They are mutable since the integrators are handled by by const-reference,
  // are initialized on-demand, and are not part of the integrator per se.

  template<CompositionMethod composition_method>
  using SRKN = SymplecticRungeKuttaNyströmIntegrator<Position,
                                                     order,
                                                     time_reversible,
                                                     evaluations,
                                                     composition_method>;
  struct NotApplicable final {};

  mutable std::conditional_t<!first_same_as_last,
                             std::unique_ptr<SRKN<BA>>,
                             NotApplicable> ba_srkn_;
  mutable std::conditional_t<first_same_as_last,
                             std::unique_ptr<SRKN<ABA>>,
                             NotApplicable> aba_srkn_;
  mutable std::conditional_t<first_same_as_last,
                             std::unique_ptr<SRKN<BAB>>,
                             NotApplicable> bab_srkn_;
};

}  // namespace internal_symplectic_runge_kutta_nyström_integrator

using internal_symplectic_runge_kutta_nyström_integrator::
    SymplecticPartitionedRungeKuttaIntegrator;

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
//   construction des tables astronomiques.
//   http://www.biodiversitylibrary.org/item/32318#page/698/mode/1up.
// - Størmer (1907), Sur les trajectoires des corpuscules électrisés dans
//   l'espace, avec application aux aurores boréales.
//   https://hal.archives-ouvertes.fr/jpa-00242574/document.
// - Verlet (1967) Computer "Experiments" on classical fluids. I.
//   Thermodynamical properties of Lennard-Jones molecules.
//   http://www.chemie.unibas.ch/~steinhauser/teaching/FS2014/MD-Simulation/Handout_12.pdf.
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/2,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/1,
                                          /*first_same_as_last=*/true> const&
NewtonDelambreStørmerVerletLeapfrog();

// Coefficients from Ruth (1983), A canonical integration technique,
// https://accelconf.web.cern.ch/accelconf/p83/PDF/PAC1983_2669.PDF.
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/3,
                                          /*time_reversible=*/false,
                                          /*evaluations=*/3,
                                          /*first_same_as_last=*/false> const&
Ruth1983();

// Coefficients from Suzuki (1990), Fractal decompositions of exponential
// operators with applications to many-body theories and Monte Carlo
// simulations.
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/5,
                                          /*first_same_as_last=*/true> const&
Suzuki1990();

// The following methods have coefficients from Yoshida (1990),
// Construction of higher order symplectic integrators
// http://sixtrack.web.cern.ch/SixTrack/doc/yoshida00.pdf.
// NOTE(egg): The coefficients were derived from equations 5.4 through 5.17
// rather than computed from the wᵢ given in tables 1 and 2.  The results were
// then cross-checked against those obtained from the tables.
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/7,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order6A();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/7,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order6B();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/7,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order6C();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8A();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8B();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8C();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8D();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
Yoshida1990Order8E();

// Coefficients from Forest and Ruth (1990),
// Fourth-order symplectic integration, equation 4.8.
// http://zwe.web.cern.ch/zwe/CAS/biblio/ruth-forest.pdf.
// This scheme was independently discovered by Candy and Rozmus (1991),
// A Symplectic Integration Algorithm for Separable Hamiltonian Functions
// (submitted earlier and published later than the Forest and Ruth paper).
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/3,
                                          /*first_same_as_last=*/true> const&
CandyRozmus1991ForestRuth1990();

// The following methods have coefficients from Robert I. McLachlan and Pau
// Atela (1992), The accuracy of symplectic integrators, table 2.
// http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/2,
                                          /*time_reversible=*/false,
                                          /*evaluations=*/2,
                                          /*first_same_as_last=*/false> const&
McLachlanAtela1992Order2Optimal();
// NOTE(egg): the coefficients given in table 2 for this integrator are
// incorrect (b1 and b2 are swapped).
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/3,
                                          /*time_reversible=*/false,
                                          /*evaluations=*/3,
                                          /*first_same_as_last=*/false> const&
McLachlanAtela1992Order3Optimal();

// The following methods have coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/2,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/2,
                                          /*first_same_as_last=*/true> const&
McLachlan1995S2();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/5,
                                          /*first_same_as_last=*/true> const&
McLachlan1995SS5();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/4,
                                          /*first_same_as_last=*/true> const&
McLachlan1995S4();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/5,
                                          /*first_same_as_last=*/true> const&
McLachlan1995S5();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/9,
                                          /*first_same_as_last=*/true> const&
McLachlan1995SS9();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/15,
                                          /*first_same_as_last=*/true> const&
McLachlan1995SS15();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/8,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/17,
                                          /*first_same_as_last=*/true> const&
McLachlan1995SS17();

// The following methods have coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/4,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/6,
                                          /*first_same_as_last=*/true> const&
BlanesMoan2002S6();
template<typename Position, typename Momentum>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          /*order=*/6,
                                          /*time_reversible=*/true,
                                          /*evaluations=*/10,
                                          /*first_same_as_last=*/true> const&
BlanesMoan2002S10();

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
