
#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "integrators/symplectic_integrator.hpp"

namespace principia {

using base::not_null;
using quantities::Quotient;

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

class SRKNIntegrator : public SymplecticIntegrator {
 protected:
  // The dimension of the time derivative.
  template<typename T>
  using Variation = Quotient<T, Time>;

 public:
  SRKNIntegrator(std::vector<double> const& a, std::vector<double> const& b);

  virtual ~SRKNIntegrator() = default;

  SRKNIntegrator() = delete;
  SRKNIntegrator(SRKNIntegrator const&) = delete;
  SRKNIntegrator(SRKNIntegrator&&) = delete;
  SRKNIntegrator& operator=(SRKNIntegrator const&) = delete;
  SRKNIntegrator& operator=(SRKNIntegrator&&) = delete;

  // The functor |compute_acceleration| computes M⁻¹ F(q, t).
  template<typename Position, typename RightHandSideComputation>
  void SolveTrivialKineticEnergyIncrement(
      RightHandSideComputation compute_acceleration,
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
  // NOTE(egg): should be std::optional.
  std::unique_ptr<FirstSameAsLast> first_same_as_last_;

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
  template<VanishingCoefficients vanishing_coefficients,
           typename Position, typename RightHandSideComputation>
  void SolveTrivialKineticEnergyIncrementOptimized(
      RightHandSideComputation compute_acceleration,
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
// Fifth order, 6 stages.  This method minimizes the error constant.  Ibidem.
SRKNIntegrator const& McLachlanAtela1992Order5Optimal();
// Sixth order, 8 stages, FSAL (synchronous momenta).
// Coefficients from Okunbor and Skeel (1994),
// Canonical Runge-Kutta-Nyström methods of orders 5 and 6,
// http://bionum.cs.purdue.edu/94OkSk.pdf.
// NOTE(egg): The coefficients were actually copied from McLachlan (1995), they
// seem to differ after a dozen significant figures or so.  Okunbor and Skeel
// remark "we did not use HYBRIDJ1 to improve the accuracy of method
// coefficients as we did in section 3.1".  We assume McLachlan's version is
// accurate.
// TODO(egg): Derive the coefficients with Mathematica.
SRKNIntegrator const& OkunborSkeel1994Order6Method13();

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_runge_kutta_nystrom_integrator_body.hpp"
