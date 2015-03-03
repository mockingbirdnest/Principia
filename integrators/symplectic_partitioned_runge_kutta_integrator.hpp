#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "integrators/symplectic_integrator.hpp"

namespace principia {

using base::not_null;

namespace integrators {

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
    double a_first;
    double a_last;
  };

  template<bool first_same_as_last,
           typename AutonomousRightHandSideComputation,
           typename RightHandSideComputation>
  void SolveOptimized(RightHandSideComputation compute_force,
                      AutonomousRightHandSideComputation compute_velocity,
                      Parameters const& parameters,
                      not_null<std::vector<SystemState>*> const solution) const;

  // Not null if the first b coefficient vanishes, in that case we can spare a
  // force computation and a velocity computation at every step.
  // NOTE(egg): should be |std::optional|.
  std::unique_ptr<FirstSameAsLast> first_same_as_last_;

  // The number of stages, or the number of stages minus one for a
  // first-same-as-last integrator.
  int stages_;

  // The position and momentum nodes.  For a first-same-as-last integrator,
  // we do not store the first b coefficient, and the last entry of |a_|
  // is the sum of the first and last a coefficients.  Use
  // |first_same_as_last_->a_first| and |first_same_as_last_->a_last| to
  // desynchronize and synchronize.
  std::vector<double> a_;
  std::vector<double> b_;

  // The weights.
  std::vector<double> c_;
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
