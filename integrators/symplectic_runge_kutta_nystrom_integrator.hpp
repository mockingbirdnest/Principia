
#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "integrators/symplectic_integrator.hpp"

namespace principia {

using base::not_null;

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

 public:
  SRKNIntegrator(std::vector<double> const& a, std::vector<double> const& b);

  virtual ~SRKNIntegrator() = default;

  SRKNIntegrator() = delete;
  SRKNIntegrator(SRKNIntegrator const&) = delete;
  SRKNIntegrator(SRKNIntegrator&&) = delete;
  SRKNIntegrator& operator=(SRKNIntegrator const&) = delete;
  SRKNIntegrator& operator=(SRKNIntegrator&&) = delete;

  // The functors |evolve_kinetic| and |evolve_potential| respectively compute
  // exp(h{·, T})(q, p) given (q, p), and exp(h{·, V})(q, p) given (q, p, t).
  template<typename Position, typename Momentum,
           typename KineticFlow,
           typename PotentialFlow>
  void SolveQuadraticKineticEnergy(
      KineticFlow evolve_kinetic,
      PotentialFlow evolve_potential,
      Parameters<Position, Momentum> const& parameters,
      not_null<Solution<Position, Momentum>*> const solution) const;

  // The functor |compute_acceleration| computes M⁻¹ F(q, t).
  template<typename Position, typename Velocity,
           typename RightHandSideComputation>
  void SolveTrivialKineticEnergyIncrement(
      RightHandSideComputation compute_acceleration,
      Parameters<Position, Velocity> const& parameters,
      not_null<Solution<Position, Velocity>*> const solution) const;

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

  VanishingCoefficients vanishing_coefficients_;
  // Null if, and only if, |vanishing_coefficients_ == kNone|.
  // If |vanishing_coefficients_ == kFirstBVanishes|, this contains the first
  // and last a coefficients.
  // If |vanishing_coefficients_ == FirstAVanishes|, this contains the first
  // and last b coefficients.
  // These are used to desynchronize and synchronize first-same-as-last
  // integrators.
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
           typename Position, typename Velocity,
           typename RightHandSideComputation>
  void SolveTrivialKineticEnergyIncrementOptimized(
      RightHandSideComputation compute_acceleration,
      Parameters<Position, Velocity> const& parameters,
      not_null<Solution<Position, Velocity>*> const solution) const;
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_runge_kutta_nystrom_integrator_body.hpp"
