// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
//  parent.
#ifndef PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
#include "integrators/ordinary_differential_equations.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_

#include <experimental/optional>
#include <list>

#include "base/status.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {

using numerics::DoublePrecision;
using numerics::FixedVector;

namespace integrators {

template<int order_>
struct AdamsMoulton {
  static int const order = order_;
  FixedVector<double, order + 1> numerators;
  double denominator;
};

//TODO(phl): Single order?
template <typename Position, int order_, int velocity_order_>
class SymmetricLinearMultistepIntegrator
    : public FixedStepSizeIntegrator<
          SpecialSecondOrderDifferentialEquation<Position>> {
  static int const half_order_ = order_ / 2 + 1;
public:
  using ODE = SpecialSecondOrderDifferentialEquation<Position>;

  SymmetricLinearMultistepIntegrator(
      serialization::FixedStepSizeIntegrator::Kind const kind,
      FixedStepSizeIntegrator<ODE> const& startup_integrator,
      AdamsMoulton<velocity_order_> const& velocity_integrator,
      FixedVector<double, half_order_> const& ɑ,
      FixedVector<double, half_order_> const& β_numerator,
      double const β_denominator);

  void Solve(Instant const& t_final,
             not_null<IntegrationInstance*> const instance) const override;

  not_null<std::unique_ptr<IntegrationInstance>> NewInstance(
    IntegrationProblem<ODE> const& problem,
    IntegrationInstance::AppendState<ODE> append_state,
    Time const& step) const;

  static int const order = order_;

 private:
  // The data for a previous step of the integration.  The |Displacement|s here
  // are really |Position|s, but we do complex computations on them and it would
  // be very inconvenient to cast these computations as barycentres.
  struct Step {
    std::vector<typename ODE::Displacement> displacements;
    std::vector<typename ODE::Acceleration> accelerations;
    DoublePrecision<Instant> time;
  };

  struct Instance : public IntegrationInstance {
    Instance(IntegrationProblem<ODE> problem,
             AppendState<ODE> append_state,
             Time step);
    ODE const equation;
    std::list<Step> previous_steps;  // At most |order_ - 1| elements.
    // During startup the following field is updated more frequently than once
    // every |step|.
    typename ODE::SystemState current_state;
    AppendState<ODE> const append_state;
    Time const step;
  };

  // Performs the startup integration, i.e., computes enough states to either
  // reach |t_final| or to reach a point where |instance.current_states| has
  // |order - 1| elements.
  void StartupSolve(Instant const& t_final,
                    Instance& instance) const;

  //TODO(phl):comment
  void VelocitySolve(int dimension,
                     Instance& instance) const;

  static void FillStepFromSystemState(ODE const& equation,
                                      typename ODE::SystemState const& state,
                                      Step& step);

  FixedStepSizeIntegrator<ODE> const& startup_integrator_;
  AdamsMoulton<velocity_order_> const& velocity_integrator_;
  FixedVector<double, half_order_> const ɑ_;
  FixedVector<double, half_order_> const β_numerator_;
  double const β_denominator_;
};

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8,
                                   /*velocity_order=*/1> const&
Quinlan1999Order8A();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8,
                                   /*velocity_order=*/4> const&
Quinlan1999Order8B();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8,
                                   /*velocity_order=*/4> const&
QuinlanTremaine1990Order8();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/10,
                                   /*velocity_order=*/4> const&
QuinlanTremaine1990Order10();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/12,
                                   /*velocity_order=*/4> const&
QuinlanTremaine1990Order12();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/14,
                                   /*velocity_order=*/4> const&
QuinlanTremaine1990Order14();

}  // namespace integrators
}  // namespace principia

#include "symmetric_linear_multistep_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
