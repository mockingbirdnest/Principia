// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
//  parent.
#ifndef PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
#include "integrators/ordinary_differential_equations.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_

#include "base/status.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {

using numerics::FixedVector;

namespace integrators {

template <typename Position, int order_>
class SymmetricLinearMultistepIntegrator
    : public FixedStepSizeIntegrator<
          SpecialSecondOrderDifferentialEquation<Position>> {
  static int const half_order_ = order_ / 2 + 1;
public:
  using ODE = SpecialSecondOrderDifferentialEquation<Position>;

  SymmetricLinearMultistepIntegrator(
      serialization::FixedStepSizeIntegrator::Kind const kind,
      FixedVector<double, half_order_> const& ɑ,
      FixedVector<double, half_order_> const& β_numerators,
      double const β_denominator);

  void Solve(Instant const& t_final,
             not_null<IntegrationInstance*> const instance) const override;

  not_null<std::unique_ptr<IntegrationInstance>> NewInstance(
    IntegrationProblem<ODE> const& problem,
    IntegrationInstance::AppendState<ODE> append_state,
    Time const& step) const;

  static int const order = order_;

 private:
  struct Instance : public IntegrationInstance {
    Instance(IntegrationProblem<ODE> problem,
             AppendState<ODE> append_state,
             Time step);
    ODE const equation;
    typename ODE::SystemState current_state;
    AppendState<ODE> const append_state;
    Time const step;
  };

  FixedVector<double, half_order_> const ɑ_;
  FixedVector<double, half_order_> β_;
};

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8> const&
Quinlan1999Order8A();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8> const&
Quinlan1999Order8B();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8> const&
QuinlanTremaine1990Order8();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/10> const&
QuinlanTremaine1990Order10();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/12> const&
QuinlanTremaine1990Order12();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/14> const&
QuinlanTremaine1990Order14();

}  // namespace integrators
}  // namespace principia

#include "symmetric_linear_multistep_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
