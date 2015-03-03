#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "integrators/symplectic_integrator.hpp"

namespace principia {

using base::not_null;

namespace integrators {

enum VanishingCoefficients {
  None,
  FirstBVanishes,
  LastAVanishes,
};

template<typename Position, typename Momentum>
class SPRKIntegrator : public SymplecticIntegrator<Position, Momentum> {
 public:
  using Coefficients = typename SymplecticIntegrator<Position,
                                                     Momentum>::Coefficients;
  using Parameters = typename SymplecticIntegrator<Position,
                                                   Momentum>::Parameters;
  using SystemState = typename SymplecticIntegrator<Position,
                                                    Momentum>::SystemState;

  SPRKIntegrator();
  ~SPRKIntegrator() override = default;

  // Coefficients from Robert I. McLachlan and Pau Atela (1992),
  // The accuracy of symplectic integrators, table 2.
  // http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf
  Coefficients const& Leapfrog() const;
  Coefficients const& Order4FirstSameAsLast() const;
  Coefficients const& Order5Optimal() const;
  // http://sixtrack.web.cern.ch/SixTrack/doc/yoshida00.pdf

  void Initialize(Coefficients const& coefficients) override;

  template<typename AutonomousRightHandSideComputation,
           typename RightHandSideComputation>
  void Solve(RightHandSideComputation compute_force,
             AutonomousRightHandSideComputation compute_velocity,
             Parameters const& parameters,
             not_null<std::vector<SystemState>*> const solution) const;

 private:

  struct FirstSameAsLast {
    double first;
    double last;
  };

  // TODO(egg): either get rid of the non-autonomous RHS, or have neither of
  // them be autonomous, this is messy.
  template<VanishingCoefficients vanishing_coefficients_,
           typename AutonomousRightHandSideComputation,
           typename RightHandSideComputation>
  void SolveOptimized(RightHandSideComputation compute_force,
                      AutonomousRightHandSideComputation compute_velocity,
                      Parameters const& parameters,
                      not_null<std::vector<SystemState>*> const solution) const;

  VanishingCoefficients vanishing_coefficients_;
  // Null if, and only if, |vanishing_coefficients_| is |None|.
  // If |vanishing_coefficients_| is |FirstBVanishes|, this contains the first
  // and last a coefficients.
  // If |vanishing_coefficients_| is |FirstAVanishes|, this contains the first
  // and last b coefficients.
  // These are used to desynchronize and synchronize first-same-as-last
  // integrators.
  std::unique_ptr<FirstSameAsLast> first_same_as_last_;

  // The number of stages, or the number of stages minus one for a
  // first-same-as-last integrator.
  int stages_;

  // The position and momentum nodes.
  // If |vanishing_coefficients_| is |FirstBVanishes|, we do not store the first
  // b coefficient, and the last entry of |a_| is the sum of the first and last
  // a coefficients.
  // If |vanishing_coefficients_| is |LastAVanishes|, we do not store the last
  // a coefficient, and the first entry of |b_| is the sum of the first and last
  // b coefficients.
  std::vector<double> a_;
  std::vector<double> b_;

  // TODO(egg): remove when we find a way to use only autonomous
  // right-hand-sides.
  // The weights.
  std::vector<double> c_;
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
