#pragma once

#include "base/macros.hpp"
#include OPTIONAL_HEADER

#include <functional>
#include <vector>

#include "base/not_null.hpp"
#include "integrators/motion_integrator.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using base::not_null;
using quantities::Time;
using quantities::Variation;

namespace integrators {

// The integrator defined in this class integrates Hamilton's equations for
// Hamiltonians of the form
//   H(q, p, t) = T(q, p) + V(q, t),
// with T quadratic in p, i.e.,
//   T(q, p) = ½ pᵀM(q)⁻¹p
// where M(q)⁻¹ is a nonsingular matrix.  The flows of T and V must be known.
// A special treatment is provided when M⁻¹(q) = M⁻¹ is independent of q and
// ∂V/∂qᵢ(q, t) is known.  In that case the differential equation can be
// reformulated as q" = M⁻¹ F(q, t), where Fᵢ = ∂V/∂qᵢ.
// TODO(egg): only the latter is implemented at this time.

// TODO(egg): a blurb on non-autonomy.

// NOTE(egg): we could implement support for time-dependent M⁻¹(q, t), see Diele
// and Marangi (2011), Explicit symplectic partitioned Runge-Kutta-Nyström
// methods for non-autonomous dynamics, but this would be a job for another
// class since it also requires a choice of a quadrature.

// TODO(egg): update the comments below to reflect the papers.
// We follow the convention of McLachlan & Atela, calling the position nodes
// aᵢ and the momentum nodes bᵢ.  The following notations appear in the
// litterature:
//   (a, b) in McLachlan and Atela, as well as Candy and Rozmus;
//   (d, c) in Ruth, Yoshida, as well as Forest and Ruth;
//   (B, b) in Sofroniou and Spaletta.
// Moreover, we follow the convention of Sofroniou and Spaletta in calling the
// weights used for the time argument of the force computation cᵢ, with
// c₁ = 0, cᵢ = cᵢ₋₁ + aᵢ₋₁ for i > 1.

class SRKNIntegrator : public MotionIntegrator {
 public:
  SRKNIntegrator(std::vector<double> const& a, std::vector<double> const& b);

  virtual ~SRKNIntegrator() = default;

  SRKNIntegrator() = delete;
  SRKNIntegrator(SRKNIntegrator const&) = delete;
  SRKNIntegrator(SRKNIntegrator&&) = delete;
  SRKNIntegrator& operator=(SRKNIntegrator const&) = delete;
  SRKNIntegrator& operator=(SRKNIntegrator&&) = delete;

  template<typename Position>
  using SRKNRightHandSideComputation =
      std::function<
          void(Time const& t,
               std::vector<Position> const&,
               not_null<std::vector<Variation<Variation<Position>>>*> const)>;

  // The functor |compute_acceleration| computes M⁻¹ F(q, t).
  template<typename Position>
  void SolveTrivialKineticEnergyIncrement(
      SRKNRightHandSideComputation<Position> compute_acceleration,
      Parameters<Position, Variation<Position>> const& parameters,
      not_null<Solution<Position, Variation<Position>>*> const solution) const;

 protected:
  enum VanishingCoefficients {
    kNone,
    kFirstBVanishes,
    kLastAVanishes,
  };

  struct FirstSameAsLast {
    double first;
    double last;
  };

  // Indicates whether some nodes vanish in a way that enables optimizations.
  VanishingCoefficients vanishing_coefficients_;
  // Null if, and only if, |vanishing_coefficients_ == kNone|.
  // If |vanishing_coefficients_ == kFirstBVanishes|, this contains the first
  // and last a coefficients.
  // If |vanishing_coefficients_ == FirstAVanishes|, this contains the first
  // and last b coefficients.
  // These are used to desynchronize and synchronize first-same-as-last
  // integrators.
  std::experimental::optional<FirstSameAsLast> first_same_as_last_;

  int stages_;

  // The position and momentum nodes.
  // If |vanishing_coefficients_ == kFirstBVanishes|, we do not store the first
  // b coefficient, and the last entry of |a_| is the sum of the first and last
  // a coefficients.
  std::vector<double> a_;
  std::vector<double> b_;

  // The weights.  Note that the first c coefficient is not stored if
  // |vanishing_coefficients_ == kFirstBVanishes|.
  std::vector<double> c_;

 private:
  template<VanishingCoefficients vanishing_coefficients, typename Position>
  void SolveTrivialKineticEnergyIncrementOptimized(
      SRKNRightHandSideComputation<Position> compute_acceleration,
      Parameters<Position, Variation<Position>> const& parameters,
      not_null<Solution<Position, Variation<Position>>*> const solution) const;
};

// Fourth order, 4 stages.  This method minimizes the error constant.
// Coefficients from Robert I. McLachlan and Pau Atela (1992),
// The accuracy of symplectic integrators, table 2.
// http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
SRKNIntegrator const& McLachlanAtela1992Order4Optimal();
// Fourth order, 5 stages, FSAL (synchronous momenta).
// Coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
SRKNIntegrator const& McLachlan1995SB3A4();
// Fourth order, 6 stages, FSAL (synchronous momenta).
// Coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
SRKNIntegrator const& McLachlan1995SB3A5();
// Fourth order, 7 stages, FSAL (synchronous positions).
// Coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
SRKNIntegrator const& BlanesMoan2002SRKN6B();
// Fifth order, 6 stages.  This method minimizes the error constant.  Ibidem.
SRKNIntegrator const& McLachlanAtela1992Order5Optimal();
// Sixth order, 8 stages, FSAL (synchronous momenta).
// Coefficients from Okunbor and Skeel (1994),
// Canonical Runge-Kutta-Nyström methods of orders 5 and 6,
// http://bionum.cs.purdue.edu/94OkSk.pdf.
// NOTE(egg): The coefficients were actually copied from McLachlan (1995), they
// seem to differ after a dozen significant figures or so.  Okunbor and Skeel
// remark "we did not use HYBRJ1 to improve the accuracy of method coefficients
// as we did in section 3.1".  We assume McLachlan's version is accurate.
// TODO(egg): Derive the coefficients with Mathematica.
SRKNIntegrator const& OkunborSkeel1994Order6Method13();
// Sixth order, 12 stages, FSAL (synchronous positions).
// Coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
SRKNIntegrator const& BlanesMoan2002SRKN11B();
// Sixth order, 15 stages, FSAL (synchronous momenta).
// Coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
SRKNIntegrator const& BlanesMoan2002SRKN14A();

}  // namespace integrators
}  // namespace principia

#include "integrators/srkn_integrator_body.hpp"
